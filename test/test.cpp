#include <iostream>

double compute(double c, double c2, double i) {
  const double T = 973.2;
  const double R = 8.314;
  const double F = 96485.33289;
  const double E_eq = 0.98636;
  const double conductivity = 6.8;
  const double I0_O2 = 0.89789;
  const double I0_H2 = 0.68601;
  const double C_O2_REF = 0.21;
  const double C_H2_REF = 1.0;
  const double C_H2O_REF = 1.0;
  const double C_R_O2 = 1;

  const double C_O_O2 = 0.21 / C_O2_REF;
  double C_R_H2 = c / C_H2_REF;
  double C_O_H2 = c2 / C_H2O_REF;
  double E = E_eq + R * T / 4. / F * log(C_R_H2);
  // const double E = E_eq + R * T / 4. / F * log(C_O_O2);
  // double eta = asinh(i / I0_O2) * R * T / F;
  double eta = 2 * R * T / F *
               log((i + sqrt(i * i + 4 * I0_H2 * I0_H2 * C_R_H2 * C_O_H2)) /
                   2. / I0_H2 / C_R_H2);
  double v = E - eta - i / conductivity;
  std::cout << "i:" << i << ",c:" << c << ",E:" << E << ",eta:" << eta
            << ",v:" << v << std::endl;
  //   std::cout << i << "\t" << v << std::endl;
  return v;
}
void scale() {
  const int WIDTH = 1020 / 1;
  const double L_ch = 0.73e-6 * 255, DX = L_ch / (WIDTH - 1);
  const double VMU_H2 = 2.0984E-5;
  const double DT = 1e-8;
  const double VMU_LAT = VMU_H2 / (DX * DX / DT);
  const double DIFF = 1.6029e-5;
  const double DIFF_LAT = DIFF / (DX * DX / DT);
  // const double VEL_IN = 0.03;
  // const double VEL_IN_LAT = 0.03 / (DX / DT);
  const double VEL_IN_LAT = 0.2;
  const double VEL_IN = VEL_IN_LAT * (DX / DT);
  std::cout << "DX:" << DX << ",DT:" << DT << ",VMU_LAT:" << VMU_LAT
            << ",DIFF_LAT:" << DIFF_LAT << ",VEL_IN:" << VEL_IN << std::endl;
  const double FTAO = 3.0 * VMU_LAT + 0.5;
  const double TS = 0.2;
  const double DTAO = 2 * DIFF_LAT / (1. - TS) + 0.5;
  std::cout << "FTAO:" << FTAO << ",DTAO:" << DTAO << std::endl;
}
int main() {
  //   compute(0.7, 1., 1.5);
  //   compute(0.9, 1., 0.5);
  // for (double i = -3; i < 3.; i += 0.1) {
  //   compute(1., 1., i);
  // }
  // compute(1, 1, 0);
  scale();
  // int a[3] = {0, -1, 3};
  // std::cout << int(uint32_t(a[1])) << std::endl;
}