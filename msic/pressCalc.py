import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import tstd
from rocketcea.cea_obj import CEA_Obj
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer, add_new_propellant


TEMP_STANDARD = 273.15
def deg_to_kelvin(temp_deg, temp_stndard=TEMP_STANDARD):
    return temp_deg + temp_stndard


def mega_to_standard(value_mega):
    return value_mega * 1000000


def calc_area_inj(r, n):
    return r ** 2 * np.pi * n


def calc_mass_flow_rate(cd, temp, p, r, n, gas_const, sigma_star):
    # TODO: create mass flow rate eq
    area_inj = calc_area_inj(r=r, n=n)
    mdot = cd * p * area_inj / np.sqrt(gas_const * temp) * sigma_star
    return mdot


def calc_mass_flow_rate_err(mfr, p, p_err, cd=1, cd_err=0):
    mfr_err = np.sqrt((mfr / p * p_err)**2 + (mfr / cd * cd_err)**2)
    return mfr_err


def calc_equivalence_ratio(mfr_f, mfr_ox, mol_mass_f, mol_mass_ox):
    fa_st = (mol_mass_f)/(3*mol_mass_ox)
    fa_exp= mfr_f/mfr_ox
    er = fa_exp/fa_st
    return er

def calc_cool_p_err(cool_p, mfr_c, mfr_c_err, mfr_tot, mfr_tot_err):
    cool_p_err = np.sqrt((cool_p/mfr_c*mfr_c_err)**2 + (cool_p/mfr_tot*mfr_tot_err)**2)
    return cool_p_err

def calc_equivalence_ratio_err(er, mfr_f, mfr_f_err, mfr_ox, mfr_ox_err):
    er_err = np.sqrt((er/mfr_f*mfr_f_err)**2 + (er/mfr_ox*mfr_ox_err)**2)
    return er_err

def calc_cstar_err(cstar, p_c, p_c_err, mfr_tot, mfr_tot_err):
    cstar_err = np.sqrt((cstar/p_c*p_c_err)**2 + (cstar/mfr_tot*mfr_tot_err)**2)
    return cstar_err

def calc_isp_err(isp, thrust, thrust_err, mfr_tot, mfr_tot_err):
    isp_err = np.sqrt((isp/thrust*thrust_err)**2 + (isp/mfr_tot*mfr_tot_err)**2)
    return isp_err


def create_df_from_rep_list(shot_name, reps_dict, df_reps=None):
    if df_reps is None:
        df_reps = pd.DataFrame()
    for rep_key, rep_value in reps_dict.items():
        df_reps.loc[shot_name, rep_key] = rep_value
    return df_reps


class RepresentativeValuesGenerater():

    def __init__(self, target_df, time_column, t_start, t_end):
        self.target_df   = target_df
        self.time_column = time_column
        self.t_start     = t_start
        self.t_end       = t_end


    def get_time_average(self, target_columnm, allow_return_tstd_err=False):
        target_values = self.target_df[(self.target_df[self.time_column] >= self.t_start) & (df1[self.time_column] <= self.t_end)][target_columnm]
        if allow_return_tstd_err:
            tstd_err = tstd(target_values, ddof=1)[0]
            return target_values.mean().iloc[0], tstd_err
        return target_values.mean().iloc[0]


class MassFlowRateParams():
    def __init__(self, target_params):
        self.cd         = target_params.get("cd",        None)
        self.cd_err     = target_params.get("cdErr",     0   )
        self.gas_const  = target_params.get("gasConst",  None)
        self.d_inj      = target_params.get("diameterM", None)
        self.n_inj      = target_params.get("numInj",    None)
        self.sigma_star = target_params.get("sigmaStar", None)
        self.mol_mass   = target_params.get("molMass",   None)


if __name__ == "__main__":

    setting = json.load(open("setting//setting.json", encoding="utf-8"))
    test_names = setting.get("testNames", [])
    df_reps_file_path = "reps//reps.csv"

    if os.path.exists(df_reps_file_path):
        df_reps = pd.read_csv(df_reps_file_path, index_col="testNames")
        print(df_reps)
    else:
        df_reps = None

    for test_name in test_names:
        # load csv
        condition_file_path = os.path.join("conditions", f"{test_name}.json")
        condition = json.load(open(condition_file_path, encoding="utf-8"))
        csv_file_path = condition.get("csvFIlePath", None)
        sensorArray = condition["sensorArray"]
        df1 = pd.read_csv(csv_file_path, names = sensorArray)
        calibration = condition.get("calibration", None)

        # Convert values to sensor calibration
        for s in sensorArray:
            if (s != "t") & (s != "trig") & (s != "T"):
                cal = calibration[sensorArray.index(s)]/10
                df1[s] = df1[s]*cal


        ox = MassFlowRateParams(condition.get("ox", None))
        f  = MassFlowRateParams(condition.get("f",  None))
        n  = MassFlowRateParams(condition.get("n",  None))

        temp_k = condition.get("tempDegC", None)
        t_start = condition.get("timeStart", None)
        t_end = condition.get("timeEnd", None)

        # calc_average_pressure
        rvg = RepresentativeValuesGenerater(
            target_df=df1,
            time_column="t",
            t_start=t_start,
            t_end=t_end
        )

        # initialize reps_dict
        reps_dict = {}

        # back_pressure Calc
        p_b_target = df1[(df1["t"] >= -2) & (df1["t"] <= 0)]["PV"]
        p_b =  p_b_target.mean()
        p_b_err = tstd(p_b_target, ddof=1)

        # mfr_calc (Ox)
        p_ox, p_ox_err = rvg.get_time_average(target_columnm=["POx"], allow_return_tstd_err=True)
        mfr_ox = calc_mass_flow_rate(
            cd=ox.cd,
            temp=deg_to_kelvin(temp_k),
            p=mega_to_standard(p_ox),
            r=ox.d_inj/2,
            n=ox.n_inj,
            gas_const=ox.gas_const,
            sigma_star=ox.sigma_star
        )
        mfr_ox_err = calc_mass_flow_rate_err(
            mfr=mfr_ox,
            p=mega_to_standard(p_ox),
            p_err=mega_to_standard(p_ox_err),
            cd=ox.cd,
            cd_err=ox.cd_err
        )
        reps_dict["POx_MPa"                       ] = p_ox
        reps_dict["POx_Error_MPa"                 ] = p_ox_err
        reps_dict["MassFlowRate_Ox_kgPerSec"      ] = mfr_ox
        reps_dict["MassFlowRate_Ox_Error_kgPerSec"] = mfr_ox_err

        # mdot_calc (F)
        p_f,  p_f_err  = rvg.get_time_average(target_columnm=["PF" ], allow_return_tstd_err=True)
        mfr_f = calc_mass_flow_rate(
            cd=f.cd,
            temp=deg_to_kelvin(temp_k),
            p=mega_to_standard(p_f),
            r=f.d_inj/2,
            n=f.n_inj,
            gas_const=f.gas_const,
            sigma_star=f.sigma_star
        )
        mfr_f_err = calc_mass_flow_rate_err(
            mfr=mfr_f,
            p=mega_to_standard(p_f),
            p_err=mega_to_standard(p_f_err),
            cd=f.cd,
            cd_err=f.cd_err
        )
        reps_dict["PF_MPa"                       ] = p_f
        reps_dict["PF_Error_MPa"                 ] = p_f_err
        reps_dict["MassFlowRate_F_kgPerSec"      ] = mfr_f
        reps_dict["MassFlowRate_F_Error_kgPerSec"] = mfr_f_err


        # mdot_calc (N)
        p_n,  p_n_err  = rvg.get_time_average(target_columnm=["PC" ], allow_return_tstd_err=True)
        mfr_n = calc_mass_flow_rate(
            cd=n.cd,
            temp=deg_to_kelvin(temp_k),
            p=mega_to_standard(p_n),
            r=n.d_inj/2,
            n=n.n_inj,
            gas_const=n.gas_const,
            sigma_star=n.sigma_star
        )
        mfr_n_err = calc_mass_flow_rate_err(
            mfr=mfr_n,
            p=mega_to_standard(p_n),
            p_err=mega_to_standard(p_n_err),
            cd=n.cd,
            cd_err=n.cd_err
        )
        filmCooling = condition.get("filmCooling", None)
        if filmCooling ==  "ON":
            reps_dict["PN_MPa"                       ] = p_n
            reps_dict["PN_Error_MPa"                 ] = p_n_err
            reps_dict["MassFlowRate_N_kgPerSec"      ] = mfr_n
            reps_dict["MassFlowRate_N_Error_kgPerSec"] = mfr_n_err
        else:
            reps_dict["PN_MPa"                       ] = p_n
            reps_dict["PN_Error_MPa"                 ] = p_n_err
            reps_dict["MassFlowRate_N_kgPerSec"      ] = 0
            reps_dict["MassFlowRate_N_Error_kgPerSec"] = 0

        # Equivalence ratio (F/OX)
        er = calc_equivalence_ratio(
            mol_mass_f=f.mol_mass,
            mol_mass_ox=ox.mol_mass,
            mfr_f=mfr_f,
            mfr_ox=mfr_ox
        )
        er_err = calc_equivalence_ratio_err(
            er=er,
            mfr_f=mfr_f,
            mfr_ox=mfr_ox,
            mfr_f_err=mfr_f_err,
            mfr_ox_err=mfr_ox_err
        )

        reps_dict["EquivalenceRatio_ul"              ] = er
        reps_dict["EquivalenceRatio_Error_ul"        ] = er_err

        filmCooling = condition.get("filmCooling", None)
        if filmCooling ==  "ON":
            #propellant
            m_prop = mfr_f + mfr_ox
            m_prop_err = mfr_f_err + mfr_ox_err

            #coolant
            m_tot = m_prop + mfr_n
            m_tot_err = m_prop_err + mfr_n_err
            cool_p = mfr_n/m_tot*100
            cool_p_err = calc_cool_p_err(
                                    cool_p=cool_p,
                                    mfr_c=mfr_n,
                                    mfr_c_err=mfr_n_err,
                                    mfr_tot= m_tot,
                                    mfr_tot_err=m_tot_err)
            reps_dict["MassFlowRate_Propellant_kgPerSec"      ] = m_prop
            reps_dict["MassFlowRate_Propellant_Error_kgPerSec"] = m_prop_err
            reps_dict["MassFlowRate_Total_kgPerSec"      ] = m_tot
            reps_dict["MassFlowRate_Total_Error_kgPerSec"] = m_tot_err
            reps_dict["Coollant_Percentage"              ] = cool_p
            reps_dict["Coollant_Percentage_Error"        ] = cool_p_err
        else:
            #propellant
            m_prop = mfr_f + mfr_ox
            m_prop_err = mfr_f_err + mfr_ox_err
            m_tot = m_prop
            m_tot_err = m_prop_err

            reps_dict["MassFlowRate_Propellant_kgPerSec"      ] = m_tot
            reps_dict["MassFlowRate_Propellant_Error_kgPerSec"] = m_tot_err
            reps_dict["MassFlowRate_Total_kgPerSec"      ] = m_tot
            reps_dict["MassFlowRate_Total_Error_kgPerSec"] = m_tot_err
            reps_dict["Coollant_Percentage"              ] = 0
            reps_dict["Coollant_Percentage_Error"        ] = 0

        # get pressure calculations:
        p_0, p_0_err = rvg.get_time_average(target_columnm=["P0"], allow_return_tstd_err=True)
        p_1, p_1_err = rvg.get_time_average(target_columnm=["P1"], allow_return_tstd_err=True)
        p_2, p_2_err = rvg.get_time_average(target_columnm=["P2"], allow_return_tstd_err=True)
        p_3, p_3_err = rvg.get_time_average(target_columnm=["P3"], allow_return_tstd_err=True)

        # write pressures
        reps_dict["P0_Pressure_MPa"                      ] = p_0
        reps_dict["P0_Error_MPa"                         ] = p_0_err
        reps_dict["P1_Pressure_MPa"                      ] = p_1
        reps_dict["P1_Error_MPa"                         ] = p_1_err
        reps_dict["P2_Pressure_MPa"                      ] = p_2
        reps_dict["P2_Error_MPa"                         ] = p_2_err
        reps_dict["P3_Pressure_MPa"                      ] = p_3
        reps_dict["P3_Error_MPa"                         ] = p_3_err

        reps_dict["Back_Pressure_MPa"                    ] = p_b
        reps_dict["Back_Pressure_Error_MPa"              ] = p_b_err

        # calculate thrust
        cal = calibration[sensorArray.index("T")]
        V_n_target = df1[(df1["t"] >= -4) & (df1["t"] <= -2)]["T"]
        V_n =  V_n_target.mean()
        V_f, V_err = rvg.get_time_average(target_columnm=["T"], allow_return_tstd_err=True)

        thrust = (cal*(V_f-V_n))
        thrust_err = (cal*(V_err))
        reps_dict["Thrust_N"                    ] = thrust
        reps_dict["Thrust_Err"                  ] = thrust_err

        #calculate ISP and C*
        At = np.pi*(0.032/2)**2
        cstar = p_1*1000000*At/m_tot
        cstar_err = calc_cstar_err(cstar = cstar,
                                    p_c = p_1,
                                    p_c_err = p_1_err,
                                    mfr_tot = m_tot,
                                    mfr_tot_err = m_tot_err)

        reps_dict["c_star_M_per_S"                      ] = cstar
        reps_dict["c_star_err_M_per_S"                      ] = cstar_err

        Isp = thrust/(m_tot*9.81)
        Isp_err = calc_isp_err(isp = Isp,
                            thrust = thrust,
                            thrust_err = thrust_err,
                            mfr_tot = m_tot,
                            mfr_tot_err = m_tot_err)
        reps_dict["ISP_s"                         ] = Isp
        reps_dict["ISP__err_s"                      ] = Isp_err

        det = condition.get("det", None)
        if det == "TRUE":
            reps_dict["Detonation"                         ] = "TRUE"
        else:
            reps_dict["Detonation"                         ] = "FALSE"

        card_str = """
        fuel C2H4  C 2.0   H 4.0
        h,cal=0.      t(k)=298.15       wt%=100.
        """
        add_new_fuel( 'C2H4', card_str )

        C = CEA_Obj(oxName="GOX", fuelName="C2H4")
        MR = C.getMRforER(ERphi = er)
        isp_theo = C.estimate_Ambient_Isp(Pc=p_1,MR=MR,Pamb=p_b,eps=p_3/p_b)[0]
        c_star_theo = C.get_Cstar(MR=MR)*0.3048

        reps_dict["c_star_exp"                      ] = cstar
        reps_dict["c_star_theo"                     ] = c_star_theo
        reps_dict["c_star_perc"                     ] = cstar/c_star_theo
        reps_dict["c_star_perc_err"                 ] = cstar_err/c_star_theo
        reps_dict["Isp_exp"                    ] = Isp
        reps_dict["Isp_theo"                   ] = isp_theo
        reps_dict["Isp_perc"                   ] = Isp/isp_theo
        reps_dict["Isp_perc_err"                    ] = Isp_err/isp_theo
        reps_dict["mass_flux"                       ] = m_tot/(np.pi*(0.032/2)**2)
        reps_dict["mass_flux_err"                   ] = m_tot_err/(np.pi*(0.032/2)**2)

        output_folder_path = os.path.join("/home/chanallen/result", test_name)
        os.makedirs(output_folder_path, exist_ok=True)
        output_file_path = os.path.join(output_folder_path, "p_t.csv")
        df1.to_csv(output_file_path, index=True)
        # raise


        # saveV_n
        df_reps = create_df_from_rep_list(test_name, reps_dict, df_reps=df_reps)
    df_reps.to_csv("reps//repsTest.csv", index_label="testNames")