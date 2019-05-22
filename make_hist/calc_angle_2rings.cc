void calc_angle_2rings(){

  TVector3 ring1(-0.110,-0.836,0.537);
  TVector3 ring2(-0.179,-0.810,0.558);
  float angle = ring1.Angle(ring2);
  cout << "angle=" << angle << endl;


}
