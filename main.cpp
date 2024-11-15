/*
 * PMC-043
 * 計算内容：一つの緯度・経度位置におけるPMCパラメタの同定を行う
 *           PMC-041は多数の再計算があるので、計算結果をメモリに入れて高速化を図る
 *
 * 実行方法：
 *
 *
 * 出力ファイル(ディレクトリは実行ディレクトリ下の data/ )：
 *  1. rayleigh_X.dat : (Xは0-2で各波長)、フィッティングする係数を求めるために計算した
 *    レイリー散乱のデータ
 *  2. obs_X.dat : フィッティングする観測データ
 *  3. optimized_pmc_X.dat : 同定されたPMCパラメタを用いて、レイリー＋PMCによる散乱を計算
 *  4. log.txt : 最適化している途中経過。誤差、
 *
 * main.cpp 緯度 経度 YYMMDD hhmm Subdirectory
 *
 *  Created on: 2024/05/27
 *      Author: ando
 */
#include <iostream>
#include <fstream>
#include <filesystem>
#include <complex>
#include <chrono>
#include <string>

#include <Msis21.h>
#include <memory_allocate.h>
#include <Mie_scattering.h>

#include "Fixed_parameter.h"
#include "Date.h"
#include "Geocoordinate.h"
#include "pmc.h"
#include "pmc_simulation.h"
#include "pmc_observation.h"
#include "msis_result.h"

std::ofstream ofs_log; /* ログ出力 */
std::ofstream ofs_param; /* 計算したパラメタの出力 */
bool logging { false };
int number_of_iteration;
bool output_log { false };

std::string process_id;

std::string obs_data_dir("/home/ando/98-Local/PMC_data/convert_data");

AndoLab::Vector3d <double> set_center_pmc(
    const double pmc_altitude, const double horizontal_shift, const double alpha);

int main(int argc, char **argv){

  wrap_msisinit_(); /* MSIS初期化 */

  /************************************************************
   * 初期化
   ************************************************************/

  Date date; /* (自作の日付クラス)解析をする日付、UT */
  double latitude, longitude; /* 解析をする(観測点の)緯度・経度 */

  AndoLab::Msis21 msis;
  std::string data_dir;

  /* 引数から、観測点の緯度・経度、(msisクラスへの)日付・時刻を取得 */
  get_arg(argc, argv, latitude, longitude, date, msis, data_dir, process_id);

  /* 連続実行のため、引数5から data下のサブディレクトリを指定する */
  std::filesystem::create_directory("data");
  std::filesystem::create_directory( data_dir );
  ofs_param.open( data_dir + "/param.txt", std::ios::out);
  ofs_log.open( data_dir + "/log.txt", std::ios::out);
  ofs_log << "# Iteration Error CenterAlt[km] Sig_z[km] Sig_r[km] R0[nm] Sig_LogNormal[-] N0[x10^6 m^-3]\n";

  /* 緯度・経度から、計算に使う座標(ひまわり方向を +x方向)への変換
   *
   * 緯度・経度・高度 0m の位置を r0 とする。地表面のtangential point
   */
  AndoLab::Vector3d <double> r0;
  convert_coordinate(latitude, longitude, r0);
  Geocoordinate Geo_r0( r0 ); /* ひまわり座標・地理座標を扱う変換 */
  std::cout << process_id << " : "
      << latitude << ", " << longitude << " ==> "
      << Geo_r0.latitude() << ", " << Geo_r0.longitude() << std::endl;

  double alpha = Geo_r0.alpha(); /* ひまわりから見た角度 */

  /************************************************************
   * 初期化（ここまで）
   ************************************************************/

  /* レイリー散乱のみの計算(こちらはフィッティングする緯度経度のみ計算) */
  double **Rayleigh = AndoLab::allocate_memory2d(Num_Lambda, N_alt, 0.0);
  double **Rayleigh2 = AndoLab::allocate_memory2d(Num_Lambda, N_alt, 0.0);

  PMC pmc[3]; /* 未設定のダミーPMCで、レイリー散乱部分を計算 (高度20-43km) */
  Msis_Result ***dummy_msis_result = new Msis_Result** [Num_Lambda];

//  std::chrono::system_clock::time_point start_time
//          = std::chrono::system_clock::now();

  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
    dummy_msis_result[j_lambda] = new Msis_Result* [N_alt];
    for(int k_alt = 0; k_alt < N_alt; k_alt++){
      dummy_msis_result[j_lambda][k_alt] = new Msis_Result [2]; /* PMC/Rayleigh interval */
    }

    std::cout << "Lambda: " << j_lambda << std::endl;
    ofs_log << "Lambda: " << j_lambda << std::endl;
    calculate_intensity(date, Lambda[j_lambda],
        Altitude_min, dAlt, N_alt,
        alpha, msis, pmc[j_lambda], Rayleigh[j_lambda],
        CALC_MSIS::RAYLEIGH_INTVL, dummy_msis_result[j_lambda]);

    /* 確認のための出力 */
    std::ofstream ofs( data_dir + "/unit_rayleigh_" + std::to_string(j_lambda) + ".dat");
    for(int i = 0; i < N_alt; i++){
      ofs << (Altitude_min + i*dAlt) * 1e-3 << " " << Rayleigh[j_lambda][i] << "\n";
    }
    ofs.close();
  }

//  std::chrono::system_clock::time_point end_time
//        = std::chrono::system_clock::now();
//  double elapsed
//   = std::chrono::duration_cast <std::chrono::milliseconds>
//   (end_time - start_time).count();
//  std::cout << "Elapsed(1) "  << elapsed * 1e-3 << " [sec]" << std::endl;
//  ofs_log << "Elapsed(1) " << elapsed * 1e-3 << " [sec]" << std::endl;
//
//  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
//    calculate_intensity(date, Lambda[j_lambda],
//        Altitude_min, dAlt, N_alt,
//        alpha, msis, pmc[j_lambda], Rayleigh[j_lambda],
//        CALC_MSIS::PMC_INTVL, dummy_msis_result[j_lambda]);
//  }
//
//  start_time
//        = std::chrono::system_clock::now();
//  elapsed
//   = std::chrono::duration_cast <std::chrono::milliseconds>
//   (start_time - end_time).count();
//  std::cout << "Elapsed(2) " << elapsed * 1e-3 << " [sec]" << std::endl;
//
//  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
//    std::cout << "Lambda: " << j_lambda << std::endl;
//    ofs_log << "Lambda: " << j_lambda << std::endl;
//    calculate_intensity(date, Lambda[j_lambda],
//        Altitude_min, dAlt, N_alt,
//        alpha, msis, pmc[j_lambda], Rayleigh[j_lambda],
//        CALC_MSIS::NO, dummy_msis_result[j_lambda]);
//
//    /* 確認のための出力 */
//    std::ofstream ofs( data_dir + "/unit_rayleigh2_" + std::to_string(j_lambda) + ".dat");
//    for(int i = 0; i < N_alt; i++){
//      ofs << (Altitude_min + i*dAlt) * 1e-3 << " " << Rayleigh[j_lambda][i] << "\n";
//    }
//    ofs.close();
//  }
//
//  end_time
//        = std::chrono::system_clock::now();
//  elapsed
//   = std::chrono::duration_cast <std::chrono::milliseconds>
//   (end_time - start_time).count();
//  std::cout << "Elapsed(3) " << elapsed * 1e-3 << " [sec]" << std::endl;
//  ofs_log << "Elapsed(3) " << elapsed * 1e-3 << " [sec]" << std::endl;
//  exit(0);

  /* 観測データの読み込み */
  /* 観測データの取得 */
  double ***Observed_data
  = AndoLab::allocate_memory3d(Num_Lambda, Num_observed_data_latitude, Num_observed_data_altitude, 0.0);
  double **Obsrvd_LatLon = AndoLab::allocate_memory2d(2, Num_observed_data_latitude, 0.0);
  get_observation_data(Observed_data, Obsrvd_LatLon, argv);

  /* 観測データのうち、この計算で求める緯度・経度のデータのインデックスを取得する */
  const int i_alpha = select_latlon( Obsrvd_LatLon, latitude, longitude );

  /* フィッティングに用いる観測データを出力する */
  for(int jLambda = 0; jLambda < Num_Lambda; jLambda++){
    std::ofstream ofs_obs( data_dir + "/obs_" + std::to_string(jLambda) + ".dat");
    for(int iAlt = 0; iAlt < Num_observed_data_altitude; iAlt++){
      ofs_obs << iAlt << " " << Observed_data[jLambda][i_alpha][iAlt] << std::endl;
    }
    ofs_obs.close();
  }

  /* フィッティング係数の導出 */
  double *SolarRayIntensity = new double [Num_Lambda];   /* 係数 = 太陽光の強さ */
  double *BackgroundIntensity = new double [Num_Lambda]; /* オフセット = 背景光の明るさ */
  ObtainFittingCoefficient(Rayleigh, i_alpha, Observed_data, SolarRayIntensity, BackgroundIntensity);

//  /* (一時的)フィッティングの様子を見るための出力 */
//  double **fit_intensity = AndoLab::allocate_memory2d(3, 100, 0.0);
//  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
//    std::ofstream ofs("data/fit_rayleigh_" + std::to_string(j_lambda) + ".dat");
//    calculate_intensity(date, Lambda[j_lambda],
//        0.0e3, dAlt, 100,
//        alpha, msis, pmc[j_lambda], fit_intensity[j_lambda]);
//
//    for(int i = 0; i < 100; i++){
////      ofs << i << " " << SolarRayIntensity[j_lambda]*fit_intensity[i] + BackgroundIntensity[j_lambda] << "\n";
//      ofs << i << " " << SolarRayIntensity[j_lambda]*fit_intensity[j_lambda][i] << "\n";
//    }
//    ofs.close();
//  }
//  /* (一時的なルーチン、ここまで) */
//
//  /* (一時的)適当な位置にPMCをおいて計算する */
//  double *temp_intensity = new double [100];
//  logging = true;
//  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
//    std::complex<double> m_ice = interpolated_m( Lambda[j_lambda] );
//    AndoLab::Vector3d <double> rc = tangential_point( 82.5246e3, alpha );
//    pmc[j_lambda].set( Lambda[j_lambda], m_ice, rc, 1.107803e3, 1500.e3, 106.811e-9, 1.11501, 6.49103e6 );
//
//    calculate_intensity(date, Lambda[j_lambda],
//        0.0e3, dAlt, 100,
//        alpha, msis, pmc[j_lambda], temp_intensity);
//    std::ofstream ofs_temp("data/temp_" + std::to_string(j_lambda) + ".dat");
//    for(int i = 0; i < 100; i++){
////      ofs_temp << i << " " << SolarRayIntensity[j_lambda]*temp_intensity[i] + BackgroundIntensity[j_lambda] << "\n";
////      ofs_temp << i << " " << SolarRayIntensity[j_lambda]*temp_intensity[i] << " "
////          << SolarRayIntensity[j_lambda]*fit_intensity[j_lambda][i] << " "
////          << SolarRayIntensity[j_lambda]*(temp_intensity[i] - fit_intensity[j_lambda][i]) << " "
////          << "\n";
//      ofs_temp << i << " " << SolarRayIntensity[j_lambda]*temp_intensity[i] + BackgroundIntensity[j_lambda] << " "
//          << SolarRayIntensity[j_lambda]*fit_intensity[j_lambda][i] << " "
//          << SolarRayIntensity[j_lambda]*(temp_intensity[i] - fit_intensity[j_lambda][i]) << " "
//          << "\n";    }
//    ofs_temp.close();
//    logging = false;
//  }
//  delete [] temp_intensity;
//  AndoLab::deallocate_memory2d(fit_intensity);
//  exit(0);
//  /* (一時的なルーチン、ここまで) */


//  /******************** デバッグのための一時的ルーチン(2024/11/5)  ********************/
//  PMC pmc_debug[3];
//  /* PMCの中心 */
//  AndoLab::Vector3d <double> rc =
//      set_center_pmc(82.5057e3, 150e3, alpha); /* PMCの中心 */
//  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
//    std::complex<double> m_ice = interpolated_m( Lambda[j_lambda] );
//    pmc_debug[j_lambda].set( Lambda[j_lambda], m_ice, rc,
//        1.10469e3, 750.0e3, 66.5536e-9, 1.39691, 22.3356e6 );
//  }
//
//  Msis_Result ***msis_debug = new Msis_Result** [Num_Lambda];
//  double **rayleigh_debug = AndoLab::allocate_memory2d(Num_Lambda, Num_PMC_Layer, 0.0);
//  double **rayleigh_debug2 = AndoLab::allocate_memory2d(Num_Lambda, Num_PMC_Layer, 0.0);
//
//  std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();
//
//  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
//    msis_debug[j_lambda] = new Msis_Result* [Num_PMC_Layer];
//    for(int k_alt = 0; k_alt < Num_PMC_Layer; k_alt++){
//      msis_debug[j_lambda][k_alt] = new Msis_Result [2];
//    }
//
//    calculate_intensity(date, Lambda[j_lambda],
//        Lower_PMC, dAlt, Num_PMC_Layer,
//        alpha, msis, pmc_debug[j_lambda], rayleigh_debug[j_lambda],
//        CALC_MSIS::PMC_INTVL, msis_debug[j_lambda]);
//    calculate_intensity(date, Lambda[j_lambda],
//        Lower_PMC, dAlt, Num_PMC_Layer,
//        alpha, msis, pmc_debug[j_lambda], rayleigh_debug2[j_lambda],
//        CALC_MSIS::RAYLEIGH_INTVL, msis_debug[j_lambda]);
//
//    std::ofstream ofs_pmc_layer(data_dir + "/debug1_" + std::to_string(j_lambda) + ".dat");
//    for(int k = 0; k < Num_PMC_Layer; k++){
//      ofs_pmc_layer << (Lower_PMC + dAlt*k) * 1e-3 << " "
//          << SolarRayIntensity[j_lambda]*rayleigh_debug[j_lambda][k] + BackgroundIntensity[j_lambda] << "\n";
//    }
//    ofs_pmc_layer.close();
//  }
//
//  std::chrono::system_clock::time_point end_time = std::chrono::system_clock::now();
//  double elapsed
//   = std::chrono::duration_cast <std::chrono::milliseconds>
//   (end_time - start_time).count();
//  std::cout << "Elapsed(A) " << elapsed * 1e-3 << " [sec]" << std::endl;
//
//  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
//    calculate_intensity(date, Lambda[j_lambda],
//        Lower_PMC, dAlt, Num_PMC_Layer,
//        alpha, msis, pmc_debug[j_lambda], rayleigh_debug[j_lambda],
//        CALC_MSIS::NO, msis_debug[j_lambda]);
//
//    std::ofstream ofs_pmc_layer(data_dir + "/debug2_" + std::to_string(j_lambda) + ".dat");
//    for(int k = 0; k < Num_PMC_Layer; k++){
//      ofs_pmc_layer << (Lower_PMC + dAlt*k) * 1e-3 << " "
//          << SolarRayIntensity[j_lambda]*rayleigh_debug[j_lambda][k] + BackgroundIntensity[j_lambda] << "\n";
//    }
//    ofs_pmc_layer.close();
//  }
//
//  start_time = std::chrono::system_clock::now();
//  elapsed
//   = std::chrono::duration_cast <std::chrono::milliseconds>
//   (start_time - end_time).count();
//  std::cout << "Elapsed(B) " << elapsed * 1e-3 << " [sec]" << std::endl;
//
//  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
//    calculate_intensity(date, Lambda[j_lambda],
//        Lower_PMC, dAlt, Num_PMC_Layer,
//        alpha, msis, pmc_debug[j_lambda], rayleigh_debug[j_lambda],
//        CALC_MSIS::YES, msis_debug[j_lambda]);
//
//    std::ofstream ofs_pmc_layer(data_dir + "/debug3_" + std::to_string(j_lambda) + ".dat");
//    for(int k = 0; k < Num_PMC_Layer; k++){
//      ofs_pmc_layer << (Lower_PMC + dAlt*k) * 1e-3 << " "
//          << SolarRayIntensity[j_lambda]*rayleigh_debug[j_lambda][k] + BackgroundIntensity[j_lambda] << "\n";
//    }
//    ofs_pmc_layer.close();
//  }
//
//  AndoLab::deallocate_memory2d(rayleigh_debug);
//  exit(0);


  /* PMC高度におけるレイリー散乱を求めておく
   * (10/16実装中)計算結果を記憶する
   */

  Msis_Result ***msis_result = new Msis_Result** [Num_Lambda];
  double **rayleigh_intensity = AndoLab::allocate_memory2d(Num_Lambda, Num_PMC_Layer, 0.0);
  double **rayleigh_intensity2 = AndoLab::allocate_memory2d(Num_Lambda, Num_PMC_Layer, 0.0);

  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){

    msis_result[j_lambda] = new Msis_Result* [Num_PMC_Layer];
    for(int k_alt = 0; k_alt < Num_PMC_Layer; k_alt++){
      msis_result[j_lambda][k_alt] = new Msis_Result [2];
    }

    calculate_intensity(date, Lambda[j_lambda],
        Lower_PMC, dAlt, Num_PMC_Layer,
        alpha, msis, pmc[j_lambda], rayleigh_intensity[j_lambda],
        CALC_MSIS::PMC_INTVL, msis_result[j_lambda]);
    calculate_intensity(date, Lambda[j_lambda],
        Lower_PMC, dAlt, Num_PMC_Layer,
        alpha, msis, pmc[j_lambda], rayleigh_intensity2[j_lambda],
        CALC_MSIS::RAYLEIGH_INTVL, msis_result[j_lambda]);

    std::ofstream ofs_pmc_layer(data_dir + "/pmc_rayleigh" + std::to_string(j_lambda) + ".dat");
    for(int k = 0; k < Num_PMC_Layer; k++){
      ofs_pmc_layer << (Lower_PMC + dAlt*k) * 1e-3 << " "
          << SolarRayIntensity[j_lambda]*rayleigh_intensity[j_lambda][k] + BackgroundIntensity[j_lambda] << "\n";
    }
    ofs_pmc_layer.close();
  }

  /* PMCパラメタの最適化をする */
  double shift_distance;
  ObtainFittingPMC(date, pmc, msis, alpha, i_alpha,
      Observed_data, rayleigh_intensity, SolarRayIntensity, BackgroundIntensity,
      shift_distance, msis_result);

//  /* 最適化されたPMCパラメタを用いて、 0-99kmを計算 */
  constexpr int Num_all_allaltitude { 100 };

  double *intensity = new double [Num_all_allaltitude];
  Msis_Result ***msis_optimized = new Msis_Result** [Num_Lambda];

//  output_log = true;

  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
    msis_optimized[j_lambda] = new Msis_Result* [ Num_all_allaltitude ];
    for(int k_alt = 0; k_alt < Num_all_allaltitude; k_alt++){
      msis_optimized[j_lambda][k_alt] = new Msis_Result [2];
    }

    calculate_intensity(date, Lambda[j_lambda],
        0.0e3, dAlt, Num_all_allaltitude,
        alpha, msis, pmc[j_lambda], intensity,
        CALC_MSIS::YES, msis_optimized[j_lambda]);

    std::ofstream ofs( data_dir + "/optimized_pmc_" + std::to_string(j_lambda) + ".dat");
    for(int i = 0; i < 100; i++){
      ofs << i << " " << SolarRayIntensity[j_lambda]*intensity[i] + BackgroundIntensity[j_lambda] << "\n";
    }
    ofs.close();
    ofs_log << "Used [" << j_lambda << "] : " << SolarRayIntensity[j_lambda] << ", " << BackgroundIntensity[j_lambda] << std::endl;
  }
  delete [] intensity;


  /* 最適化パラメタの出力 */

  /*********** tangential pointにおける高さ方向を積分したPMC数(面積あたり) **********
   *
   * shift_distance を Ds とする。tangential pointから手前方向に Ds だけ
   * 回転させたところに正規分布の中心があるので、tangential pointでは
   * exp( - Ds^2 / {2σh^2} ) になる。従って、tangential pointでの、
   * 垂直方向に積分した、単位面積あたりのPMC数 Nt は
   * Nt = N0 exp( - Ds^2 / {2σh^2} ) ∫ exp( - {r-r_0}^2 / {2σz^2} ) dr
   * = N0 exp( - Ds^2 / {2σh^2} ) √2 σz ∫ exp( - t^2 ) dt
   * = N0 exp( - Ds^2 / {2σh^2} ) √(2π) σz
   */
  const double Number_of_PMC_per_area_at_tangential_point
  = pmc[0].N0() * std::exp( - shift_distance*shift_distance / (2.0*pmc[0].sig_r()*pmc[0].sig_r() ) )
  * std::sqrt(2.0*M_PI) * pmc[0].sig_z();

  std::ofstream ofs_output( data_dir + "/output.dat" );
  ofs_output << latitude << " " << longitude << " "
      << ( pmc[0].rc().abs() - Radius_of_Earth ) * m2km << " " /* Center Altitude */
      << pmc[0].sig_z() * m2km << " " /* SD in altitude */
      << pmc[0].sig_r() * m2km << " " /* SD in altitude */
      << pmc[0].r0() * 1e9 << " " /* Particle radius */
      << pmc[0].sigma() << " "    /* SD parameter of log-normal radius distribution */
      << pmc[0].N0() * 1e-6 << " " /* Number density at the maximum (at the center) */
      << shift_distance * m2km << " " /* horizontal shift */
      << Number_of_PMC_per_area_at_tangential_point << " "
      << std::endl;
  ofs_output.close();
  ofs_log << "In main Rc(x,y,z) = "
      << pmc[0].rc().x() << ", "
      << pmc[0].rc().y() << ", "
      << pmc[0].rc().z() << std::endl;

  AndoLab::deallocate_memory2d(Rayleigh);
  AndoLab::deallocate_memory3d(Observed_data);
  AndoLab::deallocate_memory2d(Obsrvd_LatLon);
  AndoLab::deallocate_memory2d(rayleigh_intensity);
  delete [] SolarRayIntensity;
  delete [] BackgroundIntensity;


  ofs_param.close();
  ofs_log.close();

  return 0;
}


