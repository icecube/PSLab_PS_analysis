
{

  // For now, only thing this project does is load coordinate service
  // icetray library and provide public headers for interface

  cout << "Loading libpal\n";
  int fail = gSystem->Load("libpal");
  // 0 means loaded, 1 means already loaded, -1 means failed to load

  //cout << "Loading libcoordinate-service\n";
  //int fail = gSystem->Load("libcoordinate-service");
  // 0 means loaded, 1 means already loaded, -1 means failed to load
  if (fail < 0) { LOADSUCCESS = false; }
}
