/*
 * obtainfittingpmc.cpp
 *
 *  Created on: 2024/06/18
 *      Author: ando
 */
#include <iostream>
#include <fstream>
#include <vector>
//#include <chrono>
#include <string>
#include <nlopt.hpp>

#include <Msis21.h>
#include <memory_allocate.h>

#include "pmc_simulation.h"
#include "msis_result.h"


enum class OptParam : int {
  IdxAltitude = 0,
      IdxSigmaZ = 1,
      IdxSigmaH = 2,
      IdxR0 = 3,
      IdxSigma = 4,
      IdxShift = 5,
      IdxN0 = 6
};

AndoLab::Vector3d <double> rotate_alpha(AndoLab::Vector3d <double> r, const double alpha);
AndoLab::Vector3d <double> rotate_theta(AndoLab::Vector3d <double> r, const double theta);

AndoLab::Vector3d <double> set_center_pmc(
    const double pmc_altitude, const double horizontal_shift, const double alpha){

  const double R = Radius_of_Earth + pmc_altitude;
  double cos_b = R / Rgeo;
  double sin_b = std::sqrt( 1.0 - cos_b*cos_b );
  AndoLab::Vector3d <double> r(R*cos_b, 0.0, R*sin_b );
  return rotate_alpha( rotate_theta(r, horizontal_shift/R ), alpha );
}

extern std::ofstream ofs_log;

class PMC_FittingParameter{
private:
public:
  int i_alpha;
  double alpha;
  PMC *pmc;
  const Date *ptr_date;
  AndoLab::Msis21 *ptr_msis;
  const double *SolarRayIntensity;
  const double *BackgroundIntensity;
  double ***Observed_data;
  double **rayleigh_intensity;
  Msis_Result ***msis_result;
};

//constexpr double N0_base { 20.0e6 };

/* 同定に用いる評価関数。この誤差を最小化する */
double ErrPMC(const std::vector <double> &Optimized_param, std::vector <double> &grad, void *Params){

//  std::chrono::system_clock::time_point start_time
//  = std::chrono::system_clock::now();

  double **intensity = AndoLab::allocate_memory2d(Num_Lambda, Num_PMC_Layer, 0.0);

  PMC_FittingParameter *p = (PMC_FittingParameter*)Params;

  /* PMCの中心 */
  AndoLab::Vector3d <double> rc =
      set_center_pmc(Optimized_param[static_cast <int>(OptParam::IdxAltitude)],
          Optimized_param[static_cast <int>(OptParam::IdxShift)], p->alpha); /* PMCの中心 */

  double sq_err { 0.0 }; /* 二乗誤差 */

  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
    std::complex<double> m_ice = interpolated_m( Lambda[j_lambda] );
    p->pmc[j_lambda].set( Lambda[j_lambda], m_ice, rc,
        Optimized_param[static_cast <int>(OptParam::IdxSigmaZ)],
        Optimized_param[static_cast <int>(OptParam::IdxSigmaH)],
        Optimized_param[static_cast <int>(OptParam::IdxR0)],
        Optimized_param[static_cast <int>(OptParam::IdxSigma)],
        Optimized_param[static_cast <int>(OptParam::IdxN0)] );
    /* PMCの高度のみを求める */
    calculate_intensity( *( p->ptr_date ), Lambda[j_lambda],
        Lower_PMC, dAlt, Num_PMC_Layer,
        p->alpha, *(p->ptr_msis), p->pmc[j_lambda], intensity[j_lambda],
        CALC_MSIS::NO, p->msis_result[j_lambda] );

//    std::ofstream ofs("data/searching_" + std::to_string(number_of_iteration) + "_"
//        + std::to_string(j_lambda) + ".dat");
//    for(int k = 0; k < Num_PMC_Layer; k++){
//      ofs << (Lower_PMC + dAlt*k) * 1e-3 << " "
//          << p->SolarRayIntensity[j_lambda]*intensity[j_lambda][k] +
//          p->BackgroundIntensity[j_lambda] << "\n";
//    }
//    ofs.close();
  }

  /* 数密度係数を求める
   *
   */
//  double nume { 0.0 };
//  double deno { 0.0 };
//  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
//    double nume_lambda { 0.0 };
//    double deno_lambda { 0.0 };
//    for(int idx_alt = 0; idx_alt < Num_PMC_Layer; idx_alt++){
//      double I_PMC = intensity[j_lambda][idx_alt] - p->rayleigh_intensity[j_lambda][idx_alt];
//      nume_lambda += I_PMC * ( p->Observed_data[j_lambda][p->i_alpha][IdxLowerAlt + idx_alt] -
//          p->SolarRayIntensity[j_lambda] * p->rayleigh_intensity[j_lambda][idx_alt] -
//          p->BackgroundIntensity[j_lambda] );
//      deno_lambda += I_PMC * I_PMC;
//    }
//    nume += p->SolarRayIntensity[j_lambda] * nume_lambda;
//    deno += p->SolarRayIntensity[j_lambda] * p->SolarRayIntensity[j_lambda] * deno_lambda;
//  }
//  double density_coefficient { nume / deno };
//  std::cout << density_coefficient << std::endl;

//    std::ofstream ofs("data/diff_" + std::to_string(j_lambda) + ".dat");
  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
    /* 対数値の二乗誤差の計算 */
    for(int idx_alt = 0; idx_alt < Num_PMC_Layer; idx_alt++){
//      double I_PMC = intensity[j_lambda][idx_alt] - p->rayleigh_intensity[j_lambda][idx_alt];
//      double calc_result = p->SolarRayIntensity[j_lambda] *
//          ( p->rayleigh_intensity[j_lambda][idx_alt] + density_coefficient * I_PMC ) +
//          p->BackgroundIntensity[j_lambda];
      double calc_result =
          p->SolarRayIntensity[j_lambda] * intensity[j_lambda][idx_alt] +
          p->BackgroundIntensity[j_lambda];
      const int IdxLowerAlt { int( std::round( Lower_PMC / dAlt ) ) };
      double log_diffs =
          std::log( p->Observed_data[j_lambda][p->i_alpha][IdxLowerAlt + idx_alt] / calc_result );
      sq_err += log_diffs * log_diffs;
      ofs_log << "LogDiff: " << ( Lower_PMC + idx_alt*dAlt ) * m2km << " " << log_diffs
          << " " << std::log( p->Observed_data[j_lambda][p->i_alpha][IdxLowerAlt + idx_alt] )
          << " " << std::log( p->SolarRayIntensity[j_lambda] * intensity[j_lambda][idx_alt] + p->BackgroundIntensity[j_lambda] )
          << " " << p->SolarRayIntensity[j_lambda] <<
          " " << intensity[j_lambda][idx_alt] <<
          " " << p->BackgroundIntensity[j_lambda]
          << "\n";
    }
    //    ofs.close();
  }
  double SqErr = std::sqrt( sq_err );

  AndoLab::deallocate_memory2d( intensity );

//  p->pmc[0].N0( density_coefficient * N0_base ); /* 数密度係数で、N0を設定 */

  /* 最適化している様子を出力 */
  std::cout << "(" << process_id << ") "
      << number_of_iteration  << " "
      << "Err=" << SqErr << " " /* Error */
      << ( p->pmc[0].rc().abs() - Radius_of_Earth ) * m2km << "km " /* Center Altitude */
      << p->pmc[0].sig_z() * m2km << "km " /* SD in altitude */
      << p->pmc[0].sig_r() * m2km << "km " /* SD in horizontal */
      << p->pmc[0].r0() * 1e9 << "nm " /* Particle radius */
      << p->pmc[0].sigma() << " "    /* SD parameter of log-normal radius distribution */
      << p->pmc[0].N0() * 1e-6 << "cm^-3 " /* Number density at the maximum (at the center) */
      << Optimized_param[static_cast <int>(OptParam::IdxShift)] * m2km << "km " /* horizontal shift */
      << std::endl;
  ofs_log << "(" << process_id << ") "
      << number_of_iteration  << " "
      << SqErr << " " /* Error */
      << ( p->pmc[0].rc().abs() - Radius_of_Earth ) * m2km << " " /* Center Altitude */
      << p->pmc[0].sig_z() * m2km << " " /* SD in altitude */
      << p->pmc[0].sig_r() * m2km << " " /* SD in horizontal */
      << p->pmc[0].r0() * 1e9 << " " /* Particle radius */
      << p->pmc[0].sigma() << " "    /* SD parameter of log-normal radius distribution */
      << p->pmc[0].N0() * 1e-6 << " " /* Number density at the maximum (at the center) */
      << Optimized_param[static_cast <int>(OptParam::IdxShift)] * m2km << " " /* horizontal shift */
      << std::endl;
  ofs_log << "Rc(x,y,z) = "
      << p->pmc[0].rc().x() << ", "
      << p->pmc[0].rc().y() << ", "
      << p->pmc[0].rc().z() << std::endl;


//  std::chrono::system_clock::time_point end_time
//  = std::chrono::system_clock::now();
//  double elapsed
//     = std::chrono::duration_cast <std::chrono::milliseconds>
//     (end_time - start_time).count();
//  std::cout << elapsed * 1e-3 << " [sec]" << std::endl;
//  ofs_log << elapsed * 1e-3 << " [sec]" << std::endl;

  number_of_iteration++;
//  if ( number_of_iteration > Maximum_convergence ){
//    throw nlopt::forced_stop();
//  }
  return SqErr;
}


void ObtainFittingPMC(
    /*
     ********* PMCの最適化 *********
     * NLOpt を用いて PMCパラメタを同定する
     */
    Date date,
    PMC *pmc,
    AndoLab::Msis21 &msis,
    const double alpha,
    const int i_alpha,
    double ***Observed_data,
    double **rayleigh_intensity,
    const double *SolarRayIntensity,
    const double *BackgroundIntensity,
    double &shift_distance,
    Msis_Result ***msis_result){


  /*
   * アルゴリズムと、パラメタ数の設定
   * Optimized_param[0] : [m] 高度 Lower_PMC - Upper_PMC
   * Optimized_param[1] : [m] σ_z 0.1km - 5.0km
   * Optimized_param[2] : [m] σ_h (horizontal) 10km - 1500km
   * Optimized_param[3] : [m] r0 モード半径 10nm - 180nm
   * Optimized_param[4] : [-] σ 粒径分布の標準偏差 0.8-2.0
   * ×(高速化のため削除) Optimized_param[5] : [m^-3] N0 粒子数 1e5 - 200e6
   * Optimized_param[5] : [m] r_s  +y方向に対して右ねじにPMC中心位置をこの距離だけ回転(r_s = (R0 + z_pmc)*θ )
   * 0 - 750km
   *
   */
  constexpr int NUM_OPTIMIZED_PARAMETER { 7 };

  std::vector <double> Optimized_param(NUM_OPTIMIZED_PARAMETER);
  Optimized_param[static_cast <int>(OptParam::IdxAltitude)] = 82.5057e3;
  Optimized_param[static_cast <int>(OptParam::IdxSigmaZ)] = 1.10469e3;
  Optimized_param[static_cast <int>(OptParam::IdxSigmaH)] = 750.0e3;
  Optimized_param[static_cast <int>(OptParam::IdxR0)] = 66.5536e-9;
  Optimized_param[static_cast <int>(OptParam::IdxSigma)] = 1.39691;
  Optimized_param[static_cast <int>(OptParam::IdxShift)] = 150e3;
  Optimized_param[static_cast <int>(OptParam::IdxN0)] = 22.3356e6;

  /* 最大値・最小値の設定 */
  std::vector <double> LowerBound { Lower_PMC, 0.1e3, 10.0e3,    10.0e-9, 0.8, 0.0, 0.0 };
  std::vector <double> UpperBound { Upper_PMC, 5.0e3, 2000.0e3, 180.0e-9, 2.0, 1500.0e3, 1000.0e6 };

  std::vector <double> grad(NUM_OPTIMIZED_PARAMETER);

  /* わたすパラメタの設定 */
  PMC_FittingParameter Param;
  Param.alpha = alpha;
  Param.pmc = pmc;
  Param.ptr_date = &date;
  Param.ptr_msis = &msis;
  Param.SolarRayIntensity = SolarRayIntensity;
  Param.BackgroundIntensity = BackgroundIntensity;
  Param.Observed_data = Observed_data;
  Param.i_alpha = i_alpha;
  Param.rayleigh_intensity = rayleigh_intensity;
  Param.msis_result = msis_result;

//  ErrPMC( Optimized_param, grad, (void*)(&Param) ); /* お試し計算 */

  /* NLOpt による最適化 */
  nlopt::opt opt( nlopt::LN_NELDERMEAD, NUM_OPTIMIZED_PARAMETER ); /* アルゴリズムと同定パラメタ数 */
  opt.set_min_objective( ErrPMC, (void*)(&Param) ); /*　評価関数と補助パラメタ */
  opt.set_xtol_rel( 1.0e-2 ); /* Torelance */
  opt.set_lower_bounds(LowerBound);
  opt.set_upper_bounds(UpperBound);

  double minf;

  try {
    number_of_iteration = 0;
    nlopt::result result = opt.optimize( Optimized_param, minf );
  } catch (std::exception &e){
    std::cout << "Optimization of PMC parameters (" << process_id << ") failed. : " << e.what() << std::endl;
    ofs_log << "Optimization of PMC parameters (" << process_id << ") failed. : " << e.what() << std::endl;
  }

  /* pmcクラスに残っているパラメタ */
  ofs_log << "# pmc class\n"
      << number_of_iteration  << " "
      << ( pmc[0].rc().abs() - Radius_of_Earth ) * m2km << " " /* Center Altitude [km] */
      << pmc[0].sig_z() * m2km << " " /* SD in altitude [km] */
      << pmc[0].sig_r() * m2km << " " /* SD in horizontal [km] */
      << pmc[0].r0() * 1e9 << " " /* Particle radius [nm] */
      << pmc[0].sigma() << " "    /* SD parameter of log-normal radius distribution [-] */
      << pmc[0].N0() * 1e-6 << " " /* Number density at the maximum (at the center) [x10^6 m^-3] */
      << std::endl;

  /* 最適化されたパラメタ */
  ofs_log << "# opt param\n"
      << Optimized_param[static_cast <int>(OptParam::IdxAltitude)] * m2km << " "
      << Optimized_param[static_cast <int>(OptParam::IdxSigmaZ)] * m2km << " "
      << Optimized_param[static_cast <int>(OptParam::IdxSigmaH)] * m2km << " "
      << Optimized_param[static_cast <int>(OptParam::IdxR0)] * 1e9 << " "
      << Optimized_param[static_cast <int>(OptParam::IdxSigma)] << " "
      << Optimized_param[static_cast <int>(OptParam::IdxShift)] * m2km << " "
      << Optimized_param[static_cast <int>(OptParam::IdxN0)] * 1e-6 << " "
      << std::endl;

  /* 最適化された値をPMCに設定する */
//  AndoLab::Vector3d <double> rc = tangential_point( Optimized_param[0], alpha );
  AndoLab::Vector3d <double> rc =
      set_center_pmc(Optimized_param[static_cast <int>(OptParam::IdxAltitude)],
          Optimized_param[static_cast <int>(OptParam::IdxShift)], alpha);

//  const double N0 = pmc[0].N0();
  double *intensity = new double [Num_PMC_Layer]; /* デバッグ追加 */

//  PMC dummy_pmc;

  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){

    std::complex<double> m_ice = interpolated_m( Lambda[j_lambda] );
    pmc[j_lambda].set( Lambda[j_lambda], m_ice, rc,
        Optimized_param[static_cast <int>(OptParam::IdxSigmaZ)], /* σ_z */
        Optimized_param[static_cast <int>(OptParam::IdxSigmaH)], /* σ_h */
        Optimized_param[static_cast <int>(OptParam::IdxR0)], /* r0 */
        Optimized_param[static_cast <int>(OptParam::IdxSigma)], /* σ(対数正規分布) */
        Optimized_param[static_cast <int>(OptParam::IdxN0)]); /* N0 */


    calculate_intensity(date, Lambda[j_lambda],
        Lower_PMC, dAlt, Num_PMC_Layer,
//        alpha, msis, dummy_pmc, intensity,
        alpha, msis, pmc[j_lambda], intensity,
        CALC_MSIS::NO, msis_result[j_lambda]);  /* デバッグ追加 */

//    calculate_intensity(date, Lambda[j_lambda],
//        Lower_PMC, dAlt, Num_PMC_Layer,
//        alpha, msis, dummy_pmc, intensity2,
//        alpha, msis, pmc[j_lambda], intensity2,
//        CALC_MSIS::YES, msis_result[j_lambda]);  /* デバッグ追加 */
//        CALC_MSIS::YES, msis_optimized[j_lambda]);  /* デバッグ追加 */
    std::ofstream ofs( "data/optimized_pmc_insubroutine_" + std::to_string(j_lambda) + ".dat"); /* デバッグ追加 */
    for(int i = 0; i < Num_PMC_Layer; i++){ /* デバッグ追加 */
      ofs << Lower_PMC*1e-3 + i << " "
          << SolarRayIntensity[j_lambda]*intensity[i] + BackgroundIntensity[j_lambda] << "\n";
    } /* デバッグ追加 */
    ofs.close(); /* デバッグ追加 */

  }
//  delete [] msis_optimized[0];
//  delete [] msis_optimized[1];
//  delete [] msis_optimized[2];
//  delete [] msis_optimized;

  shift_distance = Optimized_param[static_cast <int>(OptParam::IdxShift)];


  delete [] intensity; /* デバッグ追加 */

}


