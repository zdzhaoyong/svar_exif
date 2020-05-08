#include <exiv2/exiv2.hpp>
#include <Svar/Svar.h>

using namespace sv;

Svar getValue(const Exiv2::Exifdatum& md){
    switch (md.typeId()) {
    case 1:
    case 3:
    case 4:
    case 6:
    case 7:
    case 8:
    case 9:
    case 14:
    case 15:
    case 16:
    case 17:
    case 18:
        if(md.count()==1)
            return (int)md.toLong();
        else
        {
            Svar ret=Svar::array();
            for(int i=0;i<md.count();i++)
                ret.push_back((int)md.toLong(i));
            return ret;
        }
    case Exiv2::signedRational:
    case Exiv2::unsignedRational:
    case Exiv2::tiffFloat:
    case Exiv2::tiffDouble:
        if(md.count()==1)
            return (double)md.toFloat();
        else
        {
            Svar ret=Svar::array();
            for(int i=0;i<md.count();i++)
                ret.push_back((double)md.toFloat(i));
            return ret;
        }
        break;
    default:
        return md.value().toString();
        break;
    };
    return Svar();
}

Svar getValue(const Exiv2::Xmpdatum& md)
{
    switch (md.typeId()) {

    case Exiv2::xmpSeq:
    {
        Svar ret=Svar::array();
        for(int i=0;i<md.count();i++)
            ret.push_back(md.toString(i));
        return ret;
    }
    default:
        return md.toString();
    };
}

Svar readMeta(std::string filepath){
    Exiv2::XmpParser::initialize();
    ::atexit(Exiv2::XmpParser::terminate);
    auto inMeta=Exiv2::ImageFactory::open(filepath);
    inMeta->readMetadata();

    Exiv2::ExifData exifIn = inMeta->exifData();
    auto xmpIn=inMeta->xmpData();

    Svar ret;
    for (Exiv2::XmpData::const_iterator md = xmpIn.begin();
         md != xmpIn.end(); ++md) {
        ret.set(md->key(),Svar({{"type",md->typeName()},{"value",getValue(*md)},{"count",(int)md->count()}}),true);
    }

    for(auto md:exifIn){
        ret.set(md.key(),Svar({{"type",md.typeName()},{"value",getValue(md)},{"count",(int)md.count()}}),true);
    }

    return ret;
}

double read_timestamp(Svar meta){
    std::string date=meta["Exif"]["Image"]["DateTime"]["value"].as<std::string>();
    int year=2020, month=0, day=0, hour=0, minute=0, second=0;
    sscanf(date.c_str(),"%d:%d:%d %d:%d:%d",&year,&month,&day,&hour,&minute,&second);
    auto TimestampFromDate=[](int year,int month,int day=0,int hour=0,int minute=0,int second=0)->double{
        tm stm;
        stm.tm_year=year-1900;
        stm.tm_mon=month-1;
        stm.tm_mday=day;
        stm.tm_hour=hour;
        stm.tm_min=minute;
        stm.tm_sec=second;
        return mktime(&stm);
    };
    return TimestampFromDate(year,month,day,hour,minute,second);
}

Svar read_gps(Svar meta){
    std::vector<double> longitude=meta["Exif"]["GPSInfo"]["GPSLongitude"]["value"].castAs<std::vector<double>>();
    double lon=longitude[0]+longitude[1]/60.+longitude[2]/3600.;

    std::vector<double> latitude=meta["Exif"]["GPSInfo"]["GPSLatitude"]["value"].castAs<std::vector<double>>();
    double lat=latitude[0]+latitude[1]/60.+latitude[2]/3600.;

    double alt=meta["Exif"]["GPSInfo"]["GPSAltitude"]["value"].castAs<double>();

    return {{"longitude",lon},{"latitude",lat},{"altitude",alt}};
}

Svar read_gpsSigma(Svar meta){
    return {{"longitude",5.},{"latitude",5.},{"altitude",10.}};
}

Svar read_attitude(const Svar& meta){
    if(meta["Xmp"]["Camera"]["Pitch"].isUndefined())
        return {{"pitch",-90.},{"yaw",0.},{"roll",0.}};

    double pitch=std::stod(meta["Xmp"]["Camera"]["Pitch"]["value"].as<std::string>())-90.;
    double yaw  =std::stod(meta["Xmp"]["Camera"]["Yaw"]["value"].as<std::string>());
    double roll =std::stod(meta["Xmp"]["Camera"]["Roll"]["value"].as<std::string>());
    return {{"pitch",pitch},{"yaw",yaw},{"roll",roll}};
}

Svar read_attitudeSigma(const Svar& meta){
    if(meta["Xmp"]["Camera"]["Pitch"].isUndefined())
        return {{"pitch",100000.},{"yaw",100000.},{"roll",100000.}};
    return {{"pitch",1.},{"yaw",10.},{"roll",1.}};
}

Svar read_camera(Svar meta)
{
    int    width=meta["Exif"]["Image"]["ImageWidth"]["value"].as<int>();
    int    height=meta["Exif"]["Image"]["ImageLength"]["value"].as<int>();
    double xres=meta["Exif"]["Photo"]["FocalPlaneXResolution"]["value"].as<double>();
    double yres=meta["Exif"]["Photo"]["FocalPlaneYResolution"]["value"].as<double>();
    double focal=std::stod(meta["Xmp"]["Camera"]["PerspectiveFocalLength"]["value"].as<std::string>());
    std::string PrincipalPoint=meta["Xmp"]["Camera"]["PrincipalPoint"]["value"].as<std::string>();
    auto    PrincipalPoints=Svar::parse_json("["+PrincipalPoint+"]");
    double cx=PrincipalPoints[0].as<double>()*xres;
    double cy=PrincipalPoints[1].as<double>()*yres;
    auto   dstr=meta["Xmp"]["Camera"]["PerspectiveDistortion"]["value"];
    std::vector<double> d;
    for(Svar s:dstr){
        std::string str=s.as<std::string>();
        d.push_back(std::stod(str));
    }

    return std::vector<double>({width,height,focal*xres,focal*yres,cx,cy,
                                d[0],d[1],d[3],d[4],d[2]});
}

#include "Matrix.h"
#include <GSLAM/core/SE3.h>

using namespace GSLAM;

Matrix3d rotation_matrix(std::vector<double> rig){
   auto deg2rad=[](double i)->double{return 3.1415925/180*i;};
   double cx = cos(deg2rad(rig[0]));
   double cy = cos(deg2rad(rig[1]));
   double cz = cos(deg2rad(rig[2]));
   double sx = sin(deg2rad(rig[0]));
   double sy = sin(deg2rad(rig[1]));
   double sz = sin(deg2rad(rig[2]));

   Matrix3d Rx ( {  1,  0,  0,
                  0, cx,-sx,
                  0, sx, cx});
   Matrix3d Ry ( { cy,  0, sy,
                  0,  1,  0,
                -sy,  0, cy});
   Matrix3d Rz ( { cz,-sz,  0,
                 sz, cz,  0,
                  0,  0,  1});
   return Rx*Ry*Rz;
}

Svar read_pose(Svar meta){
    if(meta["Xmp"]["Camera"]["RigRelatives"].isUndefined())
        return Svar();

    std::string rigss=meta["Xmp"]["Camera"]["RigRelatives"]["value"].as<std::string>();
    std::vector<double> rig=Svar::parse_json("["+rigss+"]").castAs<std::vector<double>>();
    auto m = rotation_matrix(rig);

    pi::SO3d so3;
    so3.fromMatrix(m.data());

    return {0.,0.,0.,so3.x,so3.y,so3.z,so3.w};
}

int main(int argc,char** argv){
    svar.parseMain(argc,argv);

    std::string file=svar.arg<std::string>("file","","the image file");

    std::cout<<readMeta(file)<<std::endl;
}

REGISTER_SVAR_MODULE(exiv2){
    svar["read_meta"]=readMeta;
    svar["read_timestamp"]=read_timestamp;
    svar["read_gps"] =read_gps;
    svar["read_gpsSigma"]=read_gpsSigma;
    svar["read_attitude"]=read_attitude;
    svar["read_attitudeSigma"]=read_attitudeSigma;
    svar["read_camera"]=read_camera;
    svar["read_pose"]=read_pose;
}

EXPORT_SVAR_INSTANCE
