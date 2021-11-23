#include <fstream>
#include <iomanip>
#include <iostream>

const int NX = 21, NY = 81;
int last, iter;
const double XL = 20.E-5, YL = 80.E-5;
// double x[NX + 1], y[NY + 1];
double u[NY + 2][NX + 2], v[NY + 2][NX + 2], pre[NY + 2][NX + 2];
double dx, dy, dt;
const double C = 1., CS2 = 1. / 3.;
const int FCX[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const int FCY[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
double wi[9];
double lambda[3];
double f1[NY + 2][NX + 2][9], ff1[NY + 2][NX + 2][9];
double ftao, vmu_phy, vmu_lat, preleft, preright, feq1;
int ls[NY + 2][NX + 2];
bool walls[NY + 2][NX + 2];
double delta;
double sumc_last, sumu;

void solid_structure() {
  for (int j = 0; j < NY + 2; j++) {
    for (int i = 0; i < NX + 2; i++) {
      ls[j][i] = 0;
      walls[j][i] = false;
    }
  }
  for (int i = 0; i < NX + 2; i++) {
    ls[0][i] = 1;
    ls[1][i] = 1;
    ls[NY][i] = 1;
    ls[NY + 1][i] = 1;
  }
  for (int j = 0; j < NY + 2; j++) {
    for (int i = 0; i < NX + 2; i++) {
      if (ls[j][i] == 1) {
        walls[j][i] = true;
      }
    }
  }
}
void initialization() {
  double z1, z2;
  dx = XL / (NX - 1);
  dy = dx;
  last = 500000;
  delta = 1.;
  lambda[0] = -5. / 3.;
  lambda[1] = 1. / 3.;
  lambda[2] = 1. / 12.;
  wi[0] = 4. / 9.;
  for (int i = 1; i < 5; i++) {
    wi[i] = 1. / 9.;
  }
  for (int i = 5; i < 9; i++) {
    wi[i] = 1. / 36.;
  }
  vmu_phy = 20.e-6;
  ftao = 1.;
  vmu_lat = (ftao - 0.5) / 3.;
  const double scale = vmu_phy / vmu_lat;
  dt = pow(dx, 2.) / scale;
  preleft = 1.0002;
  preright = 1.;

  for (int j = 1; j < NY + 1; j++) {
    for (int i = 1; i < NX + 1; i++) {
      pre[j][i] =
          preleft - double(i - 1) / double(NX - 1) * (preleft - preright);
      u[j][i] = 0.;
      v[j][i] = 0.;

      z2 = pow(u[j][i], 2.) + pow(v[j][i], 2.);
      for (int k = 0; k < 9; k++) {
        z1 = FCX[k] * u[j][i] + FCY[k] * v[j][i];
        const double si = wi[k] * (3. * z1 + 4.5 * pow(z1, 2.) - 1.5 * z2);
        if (k == 0) {
          feq1 = lambda[0] * pre[j][i] + si;
        } else if (k > 0 && k < 5) {
          feq1 = lambda[1] * pre[j][i] + si;
        } else if (k > 4 && k < 9) {
          feq1 = lambda[2] * pre[j][i] + si;
        }
        f1[j][i][k] = feq1;
        ff1[j][i][k] = feq1;
      }
    }
  }
  // std::cout << f1[1][2][3] << "," << ff1[1][2][3] << "," << feq1 <<
  // std::endl;
}

void collisionf() {
  double z1, z2;
  for (int j = 1; j < NY + 1; j++) {
    for (int i = 1; i < NX + 1; i++) {
      if (!walls[j][i]) {
        z2 = pow(u[j][i], 2.) + pow(v[j][i], 2.);
        for (int k = 0; k < 9; k++) {
          z1 = FCX[k] * u[j][i] + FCY[k] * v[j][i];
          const double si = wi[k] * (3. * z1 + 4.5 * pow(z1, 2.) - 1.5 * z2);
          if (k == 0) {
            feq1 = lambda[0] * pre[j][i] + si;

          } else if (k > 0 && k < 5) {
            feq1 = lambda[1] * pre[j][i] + si;
          } else if (k > 4 && k < 9) {
            feq1 = lambda[2] * pre[j][i] + si;
          }
          ff1[j][i][k] = f1[j][i][k] - 1. / ftao * (f1[j][i][k] - feq1);
        }
      }
    }
  }
}

void streamf() {
  for (int j = 1; j < NY + 1; j++) {
    for (int i = 1; i < NX + 1; i++) {
      for (int k = 0; k < 9; k++) {
        f1[j][i][k] = ff1[j - FCY[k]][i - FCX[k]][k];
      }
    }
  }
}

void boundaryf() {
  for (int j = 1; j < NY + 1; j++) {
    for (int i = 1; i < NX + 1; i++) {
      if (walls[j][i]) {
        ff1[j][i][1] = f1[j][i][3];
        ff1[j][i][3] = f1[j][i][1];
        ff1[j][i][2] = f1[j][i][4];
        ff1[j][i][4] = f1[j][i][2];
        ff1[j][i][5] = f1[j][i][7];
        ff1[j][i][7] = f1[j][i][5];
        ff1[j][i][6] = f1[j][i][8];
        ff1[j][i][8] = f1[j][i][6];
      }
    }
  }
  double z1, z2;
  for (int j = 1; j < NY + 1; j++) {
    z2 = pow(u[j][2], 2.) + pow(v[j][2], 2.);
    z1 = FCX[1] * u[j][2] + FCY[1] * v[j][2];
    feq1 = lambda[1] * pre[j][2] +
           wi[1] * (3. * z1 + 4.5 * pow(z1, 2.) - 1.5 * z2);
    f1[j][1][1] = lambda[1] * pre[j][1] +
                  wi[1] * (3. * z1 + 4.5 * pow(z1, 2.) - 1.5 * z2) +
                  f1[j][2][1] - feq1;
    z1 = FCX[5] * u[j][2] + FCY[5] * v[j][2];
    feq1 = lambda[2] * pre[j][2] +
           wi[5] * (3. * z1 + 4.5 * pow(z1, 2.) - 1.5 * z2);
    f1[j][1][5] = lambda[2] * pre[j][1] +
                  wi[5] * (3. * z1 + 4.5 * pow(z1, 2.) - 1.5 * z2) +
                  f1[j][2][5] - feq1;
    z1 = FCX[8] * u[j][2] + FCY[8] * v[j][2];
    feq1 = lambda[2] * pre[j][2] +
           wi[8] * (3. * z1 + 4.5 * pow(z1, 2.) - 1.5 * z2);
    f1[j][1][8] = lambda[2] * pre[j][1] +
                  wi[8] * (3. * z1 + 4.5 * pow(z1, 2.) - 1.5 * z2) +
                  f1[j][2][8] - feq1;

    z2 = pow(u[j][NX - 1], 2.) + pow(v[j][NX - 1], 2.);
    z1 = FCX[3] * u[j][NX - 1] + FCY[3] * v[j][NX - 1];
    feq1 = lambda[1] * pre[j][NX - 1] +
           wi[3] * (3. * z1 + 4.5 * pow(z1, 2.) - 1.5 * z2);
    f1[j][NX][3] = lambda[1] * pre[j][NX] +
                   wi[3] * (3. * z1 + 4.5 * pow(z1, 2.) - 1.5 * z2) +
                   f1[j][NX - 1][3] - feq1;
    z1 = FCX[6] * u[j][NX - 1] + FCY[6] * v[j][NX - 1];
    feq1 = lambda[2] * pre[j][NX - 1] +
           wi[6] * (3. * z1 + 4.5 * pow(z1, 2.) - 1.5 * z2);
    f1[j][NX][6] = lambda[2] * pre[j][NX] +
                   wi[6] * (3. * z1 + 4.5 * pow(z1, 2.) - 1.5 * z2) +
                   f1[j][NX - 1][6] - feq1;
    z1 = FCX[7] * u[j][NX - 1] + FCY[7] * v[j][NX - 1];
    feq1 = lambda[2] * pre[j][NX - 1] +
           wi[7] * (3. * z1 + 4.5 * pow(z1, 2.) - 1.5 * z2);
    f1[j][NX][7] = lambda[2] * pre[j][NX] +
                   wi[7] * (3. * z1 + 4.5 * pow(z1, 2.) - 1.5 * z2) +
                   f1[j][NX - 1][7] - feq1;
  }
}

void macrof() {
  for (int j = 1; j < NY + 1; j++) {
    for (int i = 1; i < NX + 1; i++) {
      if (!walls[j][i]) {
        double temppre = 0.;
        double tempu = 0.;
        double tempv = 0.;
        for (int k = 1; k < 9; k++) {
          temppre = temppre + f1[j][i][k];
          tempu = tempu + f1[j][i][k] * FCX[k];
          tempv = tempv + f1[j][i][k] * FCY[k];
        }
        u[j][i] = tempu;
        v[j][i] = tempv;
        const double temp1 = pow(u[j][i], 2.) + pow(v[j][i], 2.);
        pre[j][i] = (temppre - 2. / 3. * temp1) / (-lambda[0]);
        // std::cout << u[j][i] << "," << v[j][i] << "," << pre[j][i] <<
        // std::endl; std::endl; std::cout << temppre << "," << temp1 << "," <<
        // f1[j][i][2]
        // << std::endl;
        // std::cout << FCX[1] << "," << FCX[2] << "," << FCX[3] << std::endl;
      } else if (walls[j][i]) {
        u[j][i] = 0.;
        v[j][i] = 0.;
        pre[j][i] = 0.;
      }
    }
  }
  for (int j = 1; j < NY + 1; j++) {
    pre[j][1] = preleft;
    pre[j][NX] = preright;
  }
}

void output() {
  const double sumu_last = sumu;
  sumu = 0.;
  for (int j = 1; j < NY + 1; j++) {
    for (int i = 1; i < NX + 1; i++) {
      if (!walls[j][i]) {
        sumu = sumu + u[j][i];
        // std::cout << u[j][i] << ", ";
      }
    }
  }
  delta = abs(sumu_last - sumu) / abs(sumu);
  std::cout << iter << "," << u[NY / 2][NX - 10] << "," << delta << "," << sumu
            << std::endl;
  if (delta < 1.e-8) {
    std::ofstream outfile("velocity_pressure_cpp.dat", std::ios::out);
    outfile << "VARIABLES= X,Y,u,v,pre" << std::endl;
    outfile << "ZONE I=" << NX << ",J=" << NY << ",T=TT" << std::endl;
    for (int j = 1; j < NY + 1; j++) {
      for (int i = 1; i < NX + 1; i++) {
        outfile << i << "," << j << "," << u[j][i] << "," << v[j][i] << ","
                << pre[j][i] << std::endl;
      }
    }
    outfile.close();

    last = 0;
  }
  // last = 0;
}

int main() {
  solid_structure();
  initialization();
  for (iter = 1; iter < last + 1; iter++) {
    collisionf();
    streamf();
    boundaryf();
    macrof();
    if (iter % 1000 == 0) {
      output();
    }
  }
}