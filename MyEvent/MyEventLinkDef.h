#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class MyEvent+;
#pragma link C++ class MyTrack+;
#pragma link C++ class MyPair+;
#pragma link C++ class MyCluster+;
#pragma link C++ class std::vector<MyTrack>;
#pragma link C++ class std::vector<MyPair>;
#pragma link C++ class std::vector<MyCluster>;



#endif /* __CINT__ */