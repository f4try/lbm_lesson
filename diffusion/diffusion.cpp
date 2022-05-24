#include <math.h>

#include <fstream>
#include <iomanip>
#include <iostream>

const int NX = 1020, NY = 712;
int iter;
const double XL = 20.E-6, YL = 20.E-6;
const double DX = XL / float(NX - 1), DY = DX;
const double DIFF_PHY = 1.e-5, DTAO = 1.5, TS = 0.2,
             DIFF_LAT = (1. - TS) * (DTAO - 0.5) / 2.,
             SCALE = DIFF_PHY / DIFF_LAT, dt = pow(DX, 2.) / SCALE;
double u[NY + 1][NX + 1] = {0.}, v[NY + 1][NX + 1] = {0.},
                      con[NY + 2][NX + 2] = {0.};
const int GCX[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const int GCY[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
// const double w[9] = {4. / 9.,  1. / 9.,  1. / 9.,  1. / 9., 1. / 9.,
//                      1. / 36., 1. / 36., 1. / 36., 1. / 36.};
const double w[5] = {1. / 2., 1. / 4., 1. / 4., 1. / 4., 1. / 4.};
double g[NY + 2][NX + 2][9], gg[NY + 2][NX + 2][9];
double geq;
const double CON_IN = 1., CON_OUT = 0.;
int ls[NY + 2][NX + 2] = {0};
bool walls[NY + 2][NX + 2] = {false};
const double DELTA = 1.;
int last = 200000;
double sumc_last, sumc;

void imread(std::string filename) {
  std::ifstream imgfile(filename);
  if (!imgfile.is_open()) {
    std::cout << "can not open this file" << std::endl;
    return;
  }
  for (int y = 1; y < NY + 1; y++) {
    for (int x = 1; x < NX + 1; x++) {
      float buf;
      imgfile >> buf;
      ls[y][x] = (int)buf;
    }
  }
}
void solid_structure() {
  // imread("masks/xct_0414_left.txt");
  for (int j = 0; j < NY + 2; j++) {
    ls[j][0] = 0;
    ls[j][1] = 0;
    ls[j][2] = 0;
    ls[j][NX - 1] = 0;
    ls[j][NX] = 0;
    ls[j][NX + 1] = 0;
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
  for (int j = 1; j < NY + 1; j++) {
    for (int i = 1; i < NX + 1; i++) {
      for (int k = 0; k < 5; k++) {
        double temp_uv = GCX[k] * u[j][i] + GCY[k] * v[j][i];
        if (k == 0) {
          geq = con[j][i] * (TS + w[0] * temp_uv);
        } else {
          geq = con[j][i] * ((1. - TS) * w[k] + w[0] * temp_uv);
        }
        gg[j][i][k] = geq;
        g[j][i][k] = geq;
      }
    }
  }
}

void collisiond() {
  for (int j = 1; j < NY + 1; j++) {
    for (int i = 1; i < NX + 1; i++) {
      if (!walls[j][i]) {
        for (int k = 0; k < 5; k++) {
          const double temp_uv = GCX[k] * u[j][i] + GCY[k] * v[j][i];
          if (k == 0) {
            geq = con[j][i] * (TS + w[0] * temp_uv);
          } else {
            geq = con[j][i] * ((1. - TS) * w[k] + w[0] * temp_uv);
          }
          gg[j][i][k] = g[j][i][k] - 1. / DTAO * (g[j][i][k] - geq);
          // if (gg[j][i][k] > 0.1) {
          //   std::cout << i << "," << j << "," << k << "," << gg[j][i][k] <<
          //   ","
          //             << std::endl;
          // }
        }
      }
    }
  }
}

void streamd() {
  for (int j = 1; j < NY + 1; j++) {
    for (int i = 1; i < NX + 1; i++) {
      for (int k = 0; k < 5; k++) {
        g[j][i][k] = gg[j - GCY[k]][i - GCX[k]][k];
      }
    }
  }
}

void boundaryd() {
  for (int j = 1; j < NY + 1; j++) {
    for (int i = 1; i < NX + 1; i++) {
      if (walls[j][i]) {
        gg[j][i][1] = g[j][i][3];
        gg[j][i][3] = g[j][i][1];
        gg[j][i][2] = g[j][i][4];
        gg[j][i][4] = g[j][i][2];
        // gg[j][i][5] = g[j][i][7];
        // gg[j][i][7] = g[j][i][5];
        // gg[j][i][6] = g[j][i][8];
        // gg[j][i][8] = g[j][i][6];
      }
    }
  }
  for (int i = 1; i < NX + 1; i++) {
    gg[1][i][2] = g[1][i][4];
    gg[NY][i][4] = g[NY][i][2];
    // gg[1][i][5] = g[1][i][8];
    // gg[NY][i][8] = g[NY][i][5];
    // gg[1][i][6] = g[1][i][7];
    // gg[NY][i][7] = g[NY][i][6];
  }
  for (int j = 1; j < NY + 1; j++) {
    g[j][1][1] = CON_IN;
    g[j][NX][3] = CON_OUT;
    for (int k = 0; k < 5; k++) {
      if (k != 1) {
        g[j][1][1] -= g[j][1][k];
      }
      if (k != 3) {
        g[j][NX][3] -= g[j][NX][k];
      }
    }
  }
}

void macrod() {
  for (int j = 1; j < NY + 1; j++) {
    for (int i = 1; i < NX + 1; i++) {
      if (!walls[j][i]) {
        con[j][i] = 0.;
        for (int k = 0; k < 5; k++) {
          con[j][i] += g[j][i][k];
        }
      }
    }
  }
  for (int j = 1; j < NY + 1; j++) {
    con[j][1] = CON_IN;
    con[j][NX] = CON_OUT;
  }
}

void output() {
  const double sumc_last = sumc;
  sumc = 0.;
  for (int j = 1; j < NY + 1; j++) {
    for (int i = 1; i < NX + 1; i++) {
      // if (!walls[j][i]) {
      sumc += con[j][i];
      // }
    }
  }
  const double delta = abs(sumc_last - sumc) / abs(sumc);
  // std::cout << sumc_last << "," << sumc << "," << delta << std::endl;
  std::cout << iter << "," << iter * dt << "," << con[NY / 2][NX / 2] << ","
            << delta << "," << sumc << std::endl;
  // std::ofstream outfile("data.csv", std::ios::app);
  // outfile << iter * dt << "," << con[NY / 2][NX / 2] << std::endl;
  // outfile.close();
  // double c = 0;
  // for (int y = 2; y < NY / 5; y++) {
  //   for (int x = 2; x < NX - 2; x++) {
  //     c += con[y][x];
  //   }
  // }
  // c /= (NX - 2 - 2) * (NY / 5 - 2);
  // std::cout << c << std::endl;
  if (delta < 1.e-8) {
    std::ofstream outfile("concentration_cpp.dat", std::ios::out);
    outfile << "VARIABLES= X,Y,u,v,pre" << std::endl;
    outfile << "ZONE I=" << NX << ",J=" << NY << ",T=TT" << std::endl;
    for (int j = 2; j < NY; j++) {
      for (int i = 1; i < NX + 1; i++) {
        outfile << i << "," << j << "," << u[j][i] << "," << v[j][i] << ","
                << con[j][i] << std::endl;
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
    collisiond();
    streamd();
    boundaryd();
    macrod();
    if (iter % 100 == 0) {
      output();
    }
  }
}