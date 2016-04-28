within ;
package HySDeP

  package MHSS
    package Tanks_With_HEX

      package One_Dimensional_Models

        model Tank_system_1D

          import SI = Modelica.SIunits;

        type GravimetricDensity = Real(quantity = "Gravimetric Density", unit = "kg/kg");
        type MolarEnthalpy = Real(quantity = "MolarEnthalpy", unit = "J/mol");
        type GravimetricDensity_syst = Real(quantity = "Gravimetric Density", unit = "kg_H2/kg_syst");
        type VolumetricDensity_syst = Real(quantity = "Gravimetric Density", unit = "kg_H2/kg_syst");

          MHSS.Tanks_With_HEX.One_Dimensional_Models.MH_bed_external mH_bed_external_NEW
            annotation (Placement(transformation(extent={{-16,-2},{4,18}})));
          External_Coolant external_Coolant
            annotation (Placement(transformation(extent={{-18,-40},{2,-20}})));
          MHSS.Sources.Coolant_Source coolant_Source
            annotation (Placement(transformation(extent={{-70,-38},{-50,-18}})));
          MHSS.Tanks_With_HEX.One_Dimensional_Models.MH_bed mH_bed_NEW[N]
            annotation (Placement(transformation(extent={{-14,26},{6,46}})));

        //Ports.PressurePort p_ramp
          //  annotation (Placement(transformation(extent={{92,-10},{112,10}}),
          //      iconTransformation(extent={{92,-10},{112,10}})));
        /******************************  BOOLEAN  ***************************************************************************************************/
         inner parameter Boolean Calculated_Volume = true
            "if true then the volume is calculated from L_input and d_input";
         inner parameter Boolean VariableProperties = true
            "if true then C_MH and k_MH are interpolated from experimental values (p,Reation_Progress)";
         inner parameter Boolean Bell_Delaware_Method = true
            "If false then Kern Method";

        inner parameter Boolean Charging = true "if false then discharging";

            /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  Volumetric and Gravimetric densities definition  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        GravimetricDensity_syst Grav_Density;
        VolumetricDensity_syst Vol_Density; // it is assumed that the outer volume of the tank coincides with the inenr volume, as for a tubular tank design the vessel does not have to bear the high pressure (in the tubes) and can be made of a "thin" plastic
        SI.Mass HEX_weight;
        SI.Mass Total_H2_Mass; // total stored hydrogen mass in the tank
        SI.Mass m_s_tot; // absorbing alloy mass in whole tank

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  PARAMETERS FOR DISCRETIZATION ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        final inner parameter SI.Length L_i = L "Tank length";
        final inner parameter SI.Length dx = (delta-dx_ext)/N; // infinitesimal thickness of the external surface (so it is mre representative of the boundary where the Dirichlet condition is defined for the transferred heat to the coolant)
        final inner parameter SI.Length dx_ext = 1e-10; // infinitesimal thickness of the external surface (so it is mre representative of the boundary where the Dirichlet condition is defined for the transferred heat to the coolant)

        inner parameter Integer N = 14
            "Number of control volumes for each tube";
        // inner parameter Integer S = 4;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  BED GEOMTERY  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

         inner parameter SI.Length L = 1
            "Tank length - used if Calculated_Volume=true";
         final inner parameter SI.Length Ltube = L; // Then assumption that the tube length (HEX) equals the tank length
         inner parameter SI.Length d_tank = 0.4 "Tank diameter";
         inner parameter SI.Volume V_input = 0.180;
         inner SI.Volume V_in = (if Calculated_Volume == false then V_input else Modelica.Constants.pi*(d_tank^2/4)*L)
            "Inner volume of cylindrical tank - used if Calculated_Volume=false";

         inner SI.Length d = (if Calculated_Volume == false then (4*V_input/(Modelica.Constants.pi*L))^(1/2) else d_tank)
            "tank diameter";

        final inner parameter SI.Length ID = 2*delta
            "inner tube diameter = 2*delta";
        final inner parameter SI.Length OD = ID/0.8;
         parameter SI.Length delta=0.1e-2
            "Critical MH thickness inner tube diameter";
        inner parameter Real eps=0.6 "Bed porosity [m3_gas/m3_MH]"; // outer in MH_bed

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  COOLANT AND HEX ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        inner parameter SI.MassFlowRate mDot_c = 30 "Coolant mass flow rate";
        inner parameter SI.Temperature Tc_in = 273.15
            "Inlet coolant temperature";

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  INITIAL CONDITIONS  **************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        inner SI.Pressure p0 "Initial tank pressure";

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  PRESSURE RAMP  **************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        inner SI.Pressure p;  // =p_ramp
        inner SI.Temperature Tamb "K";
        inner SI.Pressure p_amb "Ambient pressure";
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  MH, COOLANT AND TUBE PROPERTIES  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        // Metal hydride properties stored in records [1]:
        inner replaceable parameter MHSS.Properties.MH.BaseProperties MH_properties
            constrainedby MHSS.Properties.MH.BaseProperties                     annotation (choicesAllMatching=true);
        // Assignation of values stored in the record into the variables defined in the present model
        final inner parameter GravimetricDensity w_max =  MH_properties.w_max
            "Max gravimetric density for the material [kgH2/kgTOT]";
        final inner parameter MolarEnthalpy DeltaH = MH_properties.DeltaH
            "Enthalpy of generation [J/mol]";
            final inner parameter MolarEnthalpy  DeltaH_des = MH_properties.DeltaH_des
            "Enthalpy of generation [J/mol]";
        final inner parameter SI.MolarEntropy DeltaS = MH_properties.DeltaS
            "Entropy of generation [J/(mol*K)]";
                final inner parameter SI.MolarEntropy   DeltaS_des = MH_properties.DeltaS_des
            "Entropy of generation [J/(mol*K)]";
        inner constant SI.Frequency Ca = MH_properties.Ca
            "Activation rate [1/s]";
        inner constant SI.Frequency  Ca_des = MH_properties.Ca_des
            "Activation rate [1/s]";
        inner constant MolarEnthalpy Ea = MH_properties.Ea
            "Activation Energy [J/molH2]";
            inner constant MolarEnthalpy  Ea_des = MH_properties.Ea_des
            "Activation Energy [J/molH2]";
        final inner parameter SI.Density rho_s = MH_properties.rho_s
            "Crystalline density [kg/m3]";
        final inner parameter SI.SpecificHeatCapacityAtConstantPressure cp_s0= MH_properties.cp_s0
            "Specific heat capacity [J/(kg K)]";
        final inner parameter SI.ThermalConductivity k_s = MH_properties.k_s
            "Thermal conductivity [W/(m K)]";

        inner replaceable parameter
            MHSS.Properties.Coolant.CoolantBaseProperties Coolant
            constrainedby MHSS.Properties.Coolant.CoolantBaseProperties                                         annotation (choicesAllMatching=true);

        final inner parameter SI.Density rho_c = Coolant.rho_c "kg/m3";
        final inner parameter SI.SpecificHeatCapacity cp_c = Coolant.cp_c
            "J/(kg K)";
        final inner parameter SI.ThermalConductivity k_c = Coolant.k_c "W/(m K)";
        final inner parameter SI.DynamicViscosity mu_c = Coolant.mu_c "Pa s";

        inner replaceable parameter MHSS.Properties.Tube.BaseProperties Tube_properties
            constrainedby MHSS.Properties.Tube.BaseProperties                     annotation (choicesAllMatching=true);
        //Complete_Refuelling.Properties.Tube.BaseProperties heat_exchangers(HEXCHANGER=HEXCHANGER_USED) ;
        final inner parameter SI.ThermalConductivity k_tube = Tube_properties.k_tube "W/(m K)";
        final inner parameter SI.Density rho_tube = Tube_properties.rho_tube "W/(m K)";

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  HEAT EXCHANGER  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        inner replaceable parameter
            MHSS.Heat_Transfer.Heat_Exchangers.Heat_Exchanger_Configuration HEXCHANGER_USED
            constrainedby
            MHSS.Heat_Transfer.Heat_Exchangers.Heat_Exchanger_Configuration                annotation (choicesAllMatching=true);

        /* MHSS.Heat_Transfer.Heat_Exchangers.Heat_exchangers heat_exchangers(
      HEXCHANGER=HEXCHANGER_USED)                                                                             annotation (Placement(transformation(extent={{-54,
            -132},{34,-20}})));
*/
          HySDeP.MHSS.Heat_Transfer.Heat_Exchangers.Heat_exchangers_for1D_MODEL
                                                                                 heat_exchangers_TEST_for1D(
              HEXCHANGER=HEXCHANGER_USED)
            annotation (Placement(transformation(extent={{-92,-28},{-72,-8}})));

        inner parameter Real PR = 1.25
            "Tube pitch to tube diameter ratio: Pt/OD";                             // Tube pitch to tube outer diameter ratio (Pt/OD) http://www.engineeringpage.com/technology/thermal/pitch.html
        inner parameter Modelica.SIunits.Conversions.NonSIunits.Angle_deg theta_tp = 30
            "Layout angle (30,60,45,90deg)";                                                                        // in deg
        inner parameter Integer Nt_input = 100
            "Used if NumberOfTubes_Given = true";
        inner parameter Real Nss = 1
            "Number of Sealing strips (pairs) in one baffle spacing";
        inner parameter Real B0=0.5
            "Initial value for the baffle space 0<B<1 - percentage (>0.2 (>5cm in general))";                                                  // Baffle space as percentage of the shell diameter "Ds"
        inner parameter Real Bc = 25e-2
            "Baffle cut (20-49%) with 20-25% optimum";
        inner parameter Integer Nb = 4
            "Number of baffles - input if NumberOfTubes_Given = true";
        inner parameter Real Lo_star = 1
            "1=<Lo_star<2 --> Ratio between the inlet baffle spacing and the central baffle spacing"; // Ratio between the outlet baffle spacing and the central baffle spacing
        inner parameter Real Li_star = 1
            "1=<Li_star<2 --> Ratio between the inlet baffle spacing and the central baffle spacing"; // Ratio between the inlet baffle spacing and the central baffle spacing
        inner parameter Integer N_shell = 1 "Number of shells";
        inner parameter Real Ft = 0.875
            "temperature factor (1 = countercurrent)";                             // minimum value is 0.75
        inner parameter Real a_r = 0.75 "Aspect ratio; 0<a<1 (ID/OD)";

          Adiabatic_Boundary adiabatic_Boundary
            annotation (Placement(transformation(extent={{-12,54},{8,74}})));

          Nt_distributor nt_distributor[N]
            annotation (Placement(transformation(extent={{-96,22},{-76,42}})));
          Ports.p_T p_T annotation (Placement(transformation(extent={{78,-10},{98,
                    10}}), iconTransformation(extent={{74,-10},{94,10}})));
          Ports.PressurePort p_ramp
            annotation (Placement(transformation(extent={{-102,-10},{-82,10}}),
                iconTransformation(extent={{-96,-10},{-76,10}})));
          Stored_Mass sum1
            annotation (Placement(transformation(extent={{46,-4},{66,16}})));

        equation
          m_s_tot = mH_bed_NEW[N].m_s_tube*nt_distributor[N].Nt_Port_in.Nt; //  equals the mass of absorbing alloy contained in one tube (e.g. number "N") multiplied for the number of tubes

        Grav_Density = Total_H2_Mass/(Total_H2_Mass + m_s_tot + HEX_weight); // it is assumed that the outer volume of the tank coincides with the inenr volume, as for a tubular tank design the vessel does not have to bear the high pressure (in the tubes) and can be made of a "thin" plastic
        Vol_Density = Total_H2_Mass/(V_in*1e3);
        HEX_weight = rho_tube*nt_distributor[N].Nt_Port_in.Nt*Modelica.Constants.pi*((
            OD/2)^2 - (ID/2)^2)*L;
        Total_H2_Mass = sum1.m_stored;

         p_ramp.p = p; // p_ramp coming from Pressure_ramp_and_Tamb
         p_T.T = Tamb; // ambient temperature
         p_T.p0 = p0;   // initial pressure in the tank
         p_T.p_amb = p_amb;

          connect(external_Coolant.heatTransferCoefficientPort_1D, mH_bed_external_NEW.A_Port)
            annotation (Line(
              points={{-8,-23.6},{-6,-23.6},{-6,5.4}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(coolant_Source.temperature_port, external_Coolant.C_in)
            annotation (Line(
              points={{-52.6,-28.2},{-33.3,-28.2},{-33.3,-29.9},{-15.6,-29.9}},
              color={0,0,0},
              smooth=Smooth.None));
          connect(heat_exchangers_TEST_for1D.Nt_h, mH_bed_external_NEW.Nt_h)
            annotation (Line(
              points={{-83,-12.6},{-34.5,-12.6},{-34.5,8.2},{-14.6,8.2}},
              color={0,0,0},
              smooth=Smooth.None));
         connect(mH_bed_external_NEW.B_Port, mH_bed_NEW[1].A_Port);

            for i in 1:N-1 loop
           connect(mH_bed_NEW[i].B_Port, mH_bed_NEW[i+1].A_Port);
          end for;

            for i in 1:N loop
              connect(heat_exchangers_TEST_for1D.Nt_Port, nt_distributor[i].Nt_Port_in);
              connect(nt_distributor[i].Nt_Port_out, mH_bed_NEW[i].Nt_Port);
           connect(mH_bed_NEW[i].Mass_Port, sum1.mass_Port[i]);
          end for;
        connect(mH_bed_NEW[N].B_Port, adiabatic_Boundary.C_Port);
          annotation (Diagram(coordinateSystem(extent={{-100,-100},{100,100}}),
                              graphics), Icon(coordinateSystem(extent={{-100,-100},{100,
                    100}}),                   graphics={
                  Text(
                  extent={{-26,70},{24,36}},
                  lineColor={0,0,0},
                  fillColor={170,85,255},
                  fillPattern=FillPattern.Solid,
                  textString="1D"), Bitmap(extent={{-78,-54},{78,60}}, fileName=
                     "modelica://HySDeP/Graphics/Tank1.jpg")}));
        end Tank_system_1D;

        model Adiabatic_Boundary

          Ports.HeatTransferCoefficientPort_1D C_Port
            annotation (Placement(transformation(extent={{-10,-104},{10,-84}})));

        equation
          C_Port.qDot = 0;

          annotation (Diagram(graphics), Icon(graphics={Rectangle(
                  extent={{-62,-76},{70,-32}},
                  lineColor={0,0,0},
                  fillColor={0,0,0},
                  fillPattern=FillPattern.Solid), Text(
                  extent={{-84,-8},{92,-26}},
                  lineColor={0,0,0},
                  fillColor={0,0,0},
                  fillPattern=FillPattern.Solid,
                  textString="Adiabatic boundary")}));
        end Adiabatic_Boundary;

        model MH_bed_external

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  UNITS    ***************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

          import SI = Modelica.SIunits;
          import Modelica.Constants.pi;

        type GravimetricDensity = Real(quantity = "Gravimetric Density", unit = "kg/kg");
        type MolarEnthalpy = Real(quantity = "MolarEnthalpy", unit = "J/mol");
        type MolecularWeight = Real(quantity= "MolecularWeight", unit="g/mol");
        type CoefficientOfTemperatureDerivative = Real(quantity= "CoefficientOfTemperatureDerivative", unit="J/K");
        type IsobaricThermalExpansion = Real(quantity="IsobaricThermalExpansion", unit="1/K");
        type IsothermalCompressibility = Real(quantity="IsobaricThermalExpansion", unit="1/Pa");
        type DensityDerivative = Real(quantity="DensityDerivative", unit="(kg/m3)/s");
        type DensityDerivativeAtConstantTemperature = Real(quantity="DensityDerivativeAtConstantTemperature", unit="(kg/m3)/Pa");
        type DensityDerivativeAtConstantPressure = Real(quantity="DensityDerivativeAtConstantPressure", unit="(kg/m3)/K");
        type ThermalContactResistance = Real(quantity="ThermalContactResistance", unit="(m2.K)/W");

        /************************************************  GRAPHICS  ************************************************************************************************************************************************************************************************************************************/

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  BOOLEAN  *************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        outer parameter Boolean Calculated_Volume
            "if true then the volume is calculated from L_input and d_input";

        outer parameter Boolean VariableProperties
            "if true then C_MH and k_MH are interpolated from experimental values (p,Reation_Progress)";

        outer parameter Boolean Charging "if false then discharging";

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  COOLANT AND HEX **************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
         outer parameter SI.MassFlowRate mDot_c "Coolant mass flow rate";
         outer parameter Real Ft "temperature factor (1 = countercurrent)"; // minimum value is 0.75
        SI.HeatFlowRate q_hex_ext(start=0);

        // parameter SI.CoefficientOfHeatTransfer h= 2000;
        // parameter Integer Nt = 200;
        SI.Temperature Tc_m;
        SI.ThermalResistance R_th "K/W";
        parameter ThermalContactResistance R_th_cont = 2000e-6 "[(m2.K)/W]"; //  contact resistance between tube/coolant and hydride [1]
        ThermalContactResistance R_th_tot_ext;
        SI.ThermalResistance R_th_tube "K/W";
        SI.CoefficientOfHeatTransfer U;
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  BED GEOMTERY  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

         outer parameter SI.Length L = 1
            "tank length - used if Calculated_Volume=true";
         outer parameter SI.Length d_tank = 0.4
            "tank diameter - used if Calculated_Volume=true";
         outer parameter SI.Volume V_input = 0.18
            "Inner volume of cylindrical tank - used if Calculated_Volume=false";

        SI.Length d = (if Calculated_Volume == false then (4*V_input/(Modelica.Constants.pi*L))^(1/2) else d_tank)
            "tank diameter";

        SI.Volume V_in = (if Calculated_Volume == false then V_input else Modelica.Constants.pi*(d^2/4)*L)
            "Inner volume of the cylindrical tank";

        SI.Volume V_tube_portion;

        outer parameter SI.Length L_i;
        outer parameter SI.Length dx;
        outer parameter SI.Length dx_ext; // infinitesimal thickness of the external surface (so it is mre representative of the boundary where the Dirichlet condition is defined for the transferred heat to the coolant)

         parameter Real  eps = 0.5; // porosity (Vgas/V_MH = Vgas/(Vgas+Vsolid))
        outer parameter SI.Length ID = 2*5e-3; // Inner tube inner-diameter
        final parameter SI.Length OD = ID/0.8; // Inner tube outer-diameter

        final parameter SI.Area A_in = pi*ID*L;
        final parameter SI.Area A_ext = pi*OD*L;
        SI.Area A_hex_ext;
        parameter SI.Length Local_Radius = ID/2;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  INITIAL CONDITIONS  **************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

         outer parameter SI.Temperature Tamb
            "Initial bed and hydrogen temperature";
        SI.Mass m_g0; // initial hydrogen mass stored in the tank = initial hydrogen mass in the pores (so under the assumption that at time =0 the solid is fully discharged)
        SI.Mass m_g0_full;
        SI.Mass m_des0;

        // added parameters in this test component only in order to define the variable p (i.e. pressure ramp)
        outer parameter SI.Pressure  p0 "Pa";  // Questa diventerá una inner property insieme a quella del tank (che deve essere uguale!) e verrá definita da qualche altra parte dove questo modello e Tank saranno messi dentro
        outer parameter SI.Pressure p_amb;
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  OPERATIVE VARIABLES  *************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        SI.Temp_K T_bed "Bed and hydrogen temperature";
         outer SI.Pressure p "pressure inside the tank";
        SI.Pressure p_eq "equilibrium pressure [Pa]";

        // Coefficients of the energy balance
        CoefficientOfTemperatureDerivative A
            "coefficient of temperature derivative [J/K]";
        SI.Volume B "coefficient of pressure derivative [m3]";
        SI.Heat E "coefficient of gravimetric density derivative [J]";
        SI.HeatFlowRate D "known term [W]";

        // Masses

        SI.Mass m_g_tot_tube; // total hydrogen mass stored in the tank = absorbed + gaseous in the pores [kg]
        SI.Mass m_gas_tube; // hydrogen mass stored in the pores [kg]
        SI.Mass m_abs_tube; // total hydrogen mass absorbed in the MH [kg]
        SI.MassFlowRate mDot_gas_tube;

        Real F_rp "reaction progress [kg/kg]";

        SI.HeatFlowRate q_gen_p; // Heat of compression [W]
        SI.HeatFlowRate q_gen_DeltaH; // Heat genrated during absorption reaction [W]
        SI.HeatFlowRate q_gen_tot; // Total heat generated in the tank = compession + reaction [W]

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  PROPERTIES  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        // Metal hydride properties stored in records [1]:
         outer replaceable parameter MHSS.Properties.MH.BaseProperties MH_properties
            constrainedby MHSS.Properties.MH.BaseProperties                                    annotation (choicesAllMatching=true);

        GravimetricDensity w;
        SI.Density rho_eff;  // effective density of the MH =
        SI.SpecificHeatCapacity C_MH; // effective thermal capacity = (eps*rho_gas*cp_gas + (1-eps)*rhos_s*cp_s)
        SI.ThermalConductivity k_eff; // effective thermal conductivity = (eps*k_g + (1-eps)*k_s)
         outer parameter SI.SpecificHeatCapacity cp_c;

         // Assignation of values stored in the record into the variables defined in the present model
        final parameter GravimetricDensity w_max = MH_properties.w_max
            "Max gravimetric density for the material [kgH2/kgTOT]";
        final parameter MolarEnthalpy  DeltaH = MH_properties.DeltaH
            "Enthalpy of generation [J/mol]";
            final parameter MolarEnthalpy  DeltaH_des = MH_properties.DeltaH_des
            "Enthalpy of generation [J/mol]";
        final parameter SI.MolarEntropy   DeltaS = MH_properties.DeltaS
            "Entropy of generation [J/(mol*K)]";
            final parameter SI.MolarEntropy   DeltaS_des = MH_properties.DeltaS_des
            "Entropy of generation [J/(mol*K)]";
        constant SI.Frequency  Ca = MH_properties.Ca "Activation rate [1/s]";
        constant SI.Frequency  Ca_des = MH_properties.Ca_des
            "Activation rate [1/s]";
         constant MolarEnthalpy  Ea = MH_properties.Ea
            "Activation Energy [J/molH2]";
          constant MolarEnthalpy  Ea_des = MH_properties.Ea_des
            "Activation Energy [J/molH2]";
        final parameter SI.Density   rho_s = MH_properties.rho_s
            "Crystalline density [kg/m3]";
        final parameter SI.SpecificHeatCapacityAtConstantPressure  cp_s0 = MH_properties.cp_s0
            "Specific heat capacity [J/(kg K)]";
        final parameter SI.ThermalConductivity  k_s = MH_properties.k_s
            "Thermal conductivity [W/(m K)]";
        /*

final outer parameter GravimetricDensity w_max 
    "Max gravimetric density for the material [kgH2/kgTOT]";
final outer parameter MolarEnthalpy  DeltaH "Enthalpy of generation [J/mol]";
final outer parameter SI.MolarEntropy   DeltaS 
    "Entropy of generation [J/(mol*K)]";
outer parameter SI.Frequency  Ca "Activation rate [1/s]";
 outer parameter MolarEnthalpy  Ea "Activation Energy [J/molH2]";
final outer parameter SI.Density   rho_s "Crystalline density [kg/m3]";
    final outer parameter SI.SpecificHeatCapacityAtConstantPressure  cp_s0 
    "Specific heat capacity [J/(kg K)]";
final outer parameter SI.ThermalConductivity  k_s 
    "Thermal conductivity [W/(m K)]";
*/

            // Tube properties stored in records:
        /* outer replaceable parameter Complete_Refuelling.Properties.Tube.BaseProperties Tube_properties
    constrainedby Complete_Refuelling.Properties.Tube.BaseProperties                     annotation (choicesAllMatching=true);
*/
          final outer parameter SI.ThermalConductivity  k_tube;
          // final parameter SI.ThermalConductivity  k_tube = Tube_properties.k_tube;

         constant SI.MolarEntropy R = 8.3142 "J/(mol K)";
         constant MolecularWeight MW_g = 2
            "Hydrogen molecular weight in [g/mol]";

        // hydrogen properties (CoolProp):
        replaceable package Medium =
            ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
            Modelica.Media.Interfaces.PartialMedium   annotation (choicesAllMatching=true);

         // Coolprop states and properties
         Medium.ThermodynamicState hydrogen0;
         Medium.ThermodynamicState hydrogen;
         Medium.ThermodynamicState hydrogen1;
         Medium.ThermodynamicState hydrogen2;
         Medium.ThermodynamicState hydrogen3;
         DensityDerivative drhodt;
         Real mu_g;
         SI.SpecificEnthalpy h_g;
         SI.SpecificVolume v_g;
         Real rho_dP;
         SI.SpecificHeatCapacity cp_dP;
         Real rho_dT;
         Real cv_dT;
         Real drhodP_T;
         Real drhodT_P;
         Real dcvdT_P;
         Real dcpdP_T;
         SI.Density rho_g;
         SI.Density rho_g0;
         SI.Density rho_g0_full;
         SI.SpecificInternalEnergy u_g;
         SI.SpecificHeatCapacity  cp_g;
         SI.SpecificEnthalpy h_g_in;
         SI.SpecificEntropy s_g;
         SI.SpecificHeatCapacity  cv_g;
         IsobaricThermalExpansion a_pg;
         IsothermalCompressibility k_tg;
        final parameter Real dP = 1000; // small pressure increment used in the density derivative "drhodt"
        final parameter Real dT = 0.1; //  small temperature increment used in the density derivative "drhodt"

        /*****************************************************************  LINEAR INTERPOLATION FUNCTIONS  ***************************************************************************************************/

           function DATA_interpolation
            input Real x;
            output Real y;
          protected
            Real iNew;
            Real F_data[:]={0,0.02,0.04,0.06,0.13,0.33,0.5,0.66,0.83,0.9};
            Real C_MH_data[:]={500,520,550,555,530,580,560,670,980,1050};
           algorithm
            (
           y,iNew) := Modelica.Math.Vectors.interpolate(
                 F_data,
                 C_MH_data,
                 x);
           end DATA_interpolation;

        Real x_in;
        Real y_out;

         function DATA_interpolation1
            input Real x_k_eff;
            output Real y_k_eff;
          protected
            Real iNew;
            Real p_data[:]={101300, 290000,640000,1480000,3420000,5190000,7060000,10400000,14000000,16400000,17300000,20400000,22900000,25300000};
            Real k_eff_data[:]={0.3070, 0.3077, 0.3848, 0.4606, 0.5281, 0.5572, 0.5769, 0.6030, 0.6207, 0.6534, 0.6973, 0.6672, 0.6737, 0.6677};
         algorithm
            (
           y_k_eff,iNew) := Modelica.Math.Vectors.interpolate(
                 p_data,
                 k_eff_data,
                 x_k_eff);
         end DATA_interpolation1;

        Real x_k_eff_in;
        Real y_k_eff_out;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  INITIAL EQUATIONS  ****************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

           Ports.HeatTransferCoefficientPort_1D B_Port
            annotation (Placement(transformation(extent={{-10,10},{10,30}}),
                iconTransformation(extent={{-10,22},{10,42}})));
          // Ports.Nt_h Nt_h annotation (Placement(transformation(extent={{-88,-10},{-68,10}}),
             //   iconTransformation(extent={{-88,-10},{-68,10}})));

          Ports.HeatTransferCoefficientPort_1D A_Port
            annotation (Placement(transformation(extent={{-10,-32},{10,-12}}),
                iconTransformation(extent={{-10,-36},{10,-16}})));
          Ports.Nt_h Nt_h
            annotation (Placement(transformation(extent={{-52,-10},{-32,10}}),
                iconTransformation(extent={{-96,-8},{-76,12}})));
          Ports.Mass_Port Mass_Port annotation (Placement(transformation(extent={{88,-10},
                    {108,10}}), iconTransformation(extent={{82,-8},{102,12}})));
        initial equation
          if Charging == true then
        m_g_tot_tube = m_g0;
        w = 0;
        m_abs_tube = 0;
         else
        m_g_tot_tube = m_g0_full + m_des0;
        w = 9/10*w_max;
        m_abs_tube = 9/10*m_des0;
          end if;

        T_bed = Tamb;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  EQUATIONS  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        equation

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  EQUILIBRIUM PRESSURE  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        if Charging == true then
          p_eq = p_amb * exp(DeltaH / (R * T_bed) - DeltaS / R); // when T_bed >> => p_eq >> ;  when T_bed <<  =>  p_eq <<
        else
          p_eq = p_amb * exp(DeltaH_des / (R * T_bed) - DeltaS_des / R); // when T_bed >> => p_eq >> ;  when T_bed <<  =>  p_eq <<

        end if;
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  KINETICS EQUATION  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        if Charging == true then
          if p>=p_eq then

             der(w)= Ca * exp(-Ea / (R * T_bed)) * log(p / p_eq) * (w_max - w);

           else

             der(w) = 0;

           end if;
        else
          if p>=p_eq then
            der(F_rp) = 0;
        //      der(w) = 0;

          else
               der(F_rp)= Ca_des * exp(-Ea_des / (R * T_bed)) * log(p / p_eq) * (F_rp);
          //     der(w)= Ca_des * exp(-Ea_des / (R * T_bed)) * log(p / p_eq) * (w_max - w);
          end if;
        end if;

        F_rp = w/w_max; // Reaction progress [2]

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  HYDROGEN MASS FLOW RATE  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        m_g0 = rho_g0*V_tube_portion*eps;
        m_g0_full = rho_g0_full*V_tube_portion*eps;
        m_des0 = (1 - eps) * V_tube_portion * rho_s * w_max;

        der(m_g_tot_tube) = eps * V_tube_portion * drhodt +  (1 - eps) * V_tube_portion * rho_s * der(w); // fuelled hydrogen mass flow rate (hydrogen in the pores and hydrogen  in absorbed)

        der(m_abs_tube) = (1 - eps) * V_tube_portion * rho_s * der(w);
        m_gas_tube = eps * V_tube_portion * rho_g;
        mDot_gas_tube = eps * V_tube_portion * drhodt;

        /*
m_g0 = rho_g0*V_tube_portion*eps;
der(m_g_tot_tube) = eps * V_tube_portion * drhodt +  (1 - eps) * V_tube_portion * rho_s * der(w); // fuelled hydrogen mass flow rate (hydrogen in the pores and hydrogen  in absorbed)

der(m_gas_tube) = eps * V_tube_portion * drhodt;
der(m_abs_tube) = (1 - eps) * V_tube_portion * rho_s * der(w);
*/
        /***************************************************************************************************************************************************************
******************************************************************  SOLID BED TEMPERATURE  *********************************************************************
****************************************************************************************************************************************************************/

         A = ((1-eps)*rho_s*w*cp_g + rho_eff*C_MH)*V_tube_portion;  // Coefficient of Temperature derivative - Considers the effective heat capacity C_MH = eps*rho_g*cp_g + (1-eps)*rho_s*cp_s = 500 J/(kg m3)

         B = ((1-eps) + ((1-eps)*rho_s/rho_g*w+eps)*(1 - a_pg*T_bed) - eps)*V_tube_portion; // Coefficient of pressure derivative

        if Charging == true then
        D = (der(m_g_tot_tube)*(h_g_in - h_g) + A_Port.qDot + B_Port.qDot); // known term (energy in, and its time variation in the tank, and exchanged heat)
         // D = (der(m_g_tot_tube)*(h_g_in-h_g) + A_Port.qDot + B_Port.qDot); // known term (energy in, and its time variation in the tank, and exchanged heat)

         E = ((1-eps)*rho_s * DeltaH/(MW_g/1000))*V_tube_portion;  // Coefficient of gravimetric-density derivative

        else
          D = (der(m_g_tot_tube)*(h_g_in - h_g) + A_Port.qDot + B_Port.qDot);    // known term (energy in, and its time variation in the tank, and exchanged heat)
          E = ((1-eps)*rho_s * DeltaH_des/(MW_g/1000))*V_tube_portion;  // Coefficient of gravimetric-density derivative

        end if;

         V_tube_portion = pi*((ID/2)^2 - (ID/2-dx_ext)^2)*L_i*Nt_h.Nt;

        q_hex_ext = U*A_hex_ext*(T_bed - Tc_m); // heat exchanged for one tube only
        A_hex_ext = pi*ID*L_i*Nt_h.Nt;
         U = 1/R_th_tot_ext;
         R_th_tot_ext = (1/(1/Nt_h.h + R_th_cont + (R_th*A_hex_ext) + (R_th_tube*A_hex_ext)))^(-1); // [m2 W/K] => Rth*A_hex_ext = dx_ext/2/k_eff [m2 W/K]
          R_th_tube = OD/2*log(OD/ID)/(A_hex_ext*k_tube); // [K/W] thermal resistance of the tube wall

        R_th = dx_ext/2/(A_hex_ext*k_eff); // [K/W] thermal resistance with respect to half thickness MUST BE DIVIDED BY 2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        A_Port.T = Tc_m;
        A_Port.qDot = -q_hex_ext;
        A_Port.Counter = 0;

        B_Port.Counter = A_Port.Counter +1;
        B_Port.qDot = (B_Port.T-T_bed)/(R_th); // the effect of the number of tubes is already included in the definition of A_hex

        Mass_Port.m_tot = m_g_tot_tube;

        der(T_bed)= (D - B*der(p) - E *der(w))/A; // Energy balance

        /***************************************************************************************************************************************************************
  ***************************************************************  MH properties  ******************************************************************************
  **************************************************************************************************************************************************************/
        rho_eff = (1-eps)*rho_s + eps * rho_g; // effective density

        x_in = F_rp;
        y_out = DATA_interpolation(x=x_in);

        x_k_eff_in = p;
        y_k_eff_out = DATA_interpolation1(x_k_eff=x_k_eff_in);

        if VariableProperties == true then

        C_MH = y_out;
        k_eff = y_k_eff_out;

        else
        k_eff = MH_properties.k_s;
          C_MH = 500;
        end if;

        /***************************************************************************************************************************************************************
  ***************************************************************  Heat in the tank  ******************************************************************************
  **************************************************************************************************************************************************************/

        q_gen_p = eps*V_tube_portion*der(p); // Heat of compression [1]
        q_gen_DeltaH = -E*der(w); // Heat of reaction
        q_gen_tot = q_gen_p + q_gen_DeltaH;

        /***************************************************************************************************************************************************************
  *******************************************  Thermodynamic states and properties calculations (COOLPROP)  ****************************************************
  **************************************************************************************************************************************************************/

         drhodt = drhodP_T * der(p) + drhodT_P * der(T_bed); // derivative of the hydrogen density

          v_g=1/rho_g;

         hydrogen0=Medium.setState_pT(p=p0, T=Tamb);
         hydrogen=Medium.setState_pT(p=p, T=T_bed);
         hydrogen1=Medium.setState_pT(p=p, T=Tamb);
         hydrogen2=Medium.setState_pT(p=p+dP, T=T_bed);
         hydrogen3=Medium.setState_pT(p=p, T=T_bed+dT);

         h_g_in=Medium.specificEnthalpy(hydrogen1);
         h_g=Medium.specificEnthalpy(hydrogen) "Dynamic enthalpy";
         cp_g=Medium.specificHeatCapacityCp(hydrogen);
         u_g=Medium.specificInternalEnergy(hydrogen);
          rho_g0=Medium.density(hydrogen0);
         rho_g=Medium.density(hydrogen);
           rho_g0_full=Medium.density(hydrogen0);
         s_g=Medium.specificEntropy(hydrogen);
         cv_g=Medium.specificHeatCapacityCv(hydrogen);
         mu_g=Medium.dynamicViscosity(hydrogen);

         rho_dP =Medium.density(hydrogen2);
         cp_dP= Medium.specificHeatCapacityCp(hydrogen2);

         drhodP_T = (rho_dP - rho_g) / dP;
         drhodT_P = (rho_dT - rho_g) / dT;

         a_pg = rho_g*(1/rho_dT - 1/rho_g) / dT;
         k_tg = - rho_g*(1/rho_dP - 1/rho_g)/dP;

         rho_dT = Medium.density(hydrogen3);
         cv_dT = Medium.specificHeatCapacityCv(hydrogen3);

         dcpdP_T = (cp_dP - cp_g) / dP;
         dcvdT_P = (cv_dT - cv_g) / dT;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************ REFERENCES ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /*
[1] M. Visaria, I. Mudawar, T. Pourpoint, and S. Kumar, “Study of heat transfer and kinetics parameters influencing 
    the design of heat exchangers for hydrogen storage in high-pressure metal hydrides,” Int. J. Heat Mass Transf., vol. 
    53, no. 9–10, pp. 2229–2239, Apr. 2010.
[2] T. L. Pourpoint, V. Velagapudi, I. Mudawar, Y. Zheng, and T. S. Fisher, “Active cooling of a metal hydride 
    system for hydrogen storage,” Int. J. Heat Mass Transf., vol. 53, no. 7–8, pp. 1326–1332, Mar. 2010.
[3]
[4]

*/

         annotation (Diagram(graphics), Icon(graphics={Bitmap(extent={{-72,-52},
                      {82,56}},  fileName=
                      "modelica://HySDeP/../../Figures for components/MH_bed.png"),                     Text(
                  extent={{-98,-18},{-16,-50}},
                  lineColor={0,0,255},
                  fillColor={255,170,85},
                  fillPattern=FillPattern.Solid,
                  textString="Ext")}));
        end MH_bed_external;

        model MH_bed

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  UNITS    ***************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

          import SI = Modelica.SIunits;
          import Modelica.Constants.pi;

        type GravimetricDensity = Real(quantity = "Gravimetric Density", unit = "kg/kg");
        type MolarEnthalpy = Real(quantity = "MolarEnthalpy", unit = "J/mol");
        type MolecularWeight = Real(quantity= "MolecularWeight", unit="g/mol");
        type CoefficientOfTemperatureDerivative = Real(quantity= "CoefficientOfTemperatureDerivative", unit="J/K");
        type IsobaricThermalExpansion = Real(quantity="IsobaricThermalExpansion", unit="1/K");
        type IsothermalCompressibility = Real(quantity="IsobaricThermalExpansion", unit="1/Pa");
        type DensityDerivative = Real(quantity="DensityDerivative", unit="(kg/m3)/s");
        type DensityDerivativeAtConstantTemperature = Real(quantity="DensityDerivativeAtConstantTemperature", unit="(kg/m3)/Pa");
        type DensityDerivativeAtConstantPressure = Real(quantity="DensityDerivativeAtConstantPressure", unit="(kg/m3)/K");
        /************************************************  GRAPHICS  ************************************************************************************************************************************************************************************************************************************/

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  BOOLEAN  *************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

         outer parameter Boolean Calculated_Volume
            "if true then the volume is calculated from L_input and d_input";

         outer parameter Boolean VariableProperties
            "if true then C_MH and k_MH are interpolated from experimental values (p,Reation_Progress)";

        outer parameter Boolean Charging "if false then discharging";

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  COOLANT AND HEX **************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
         outer parameter SI.MassFlowRate mDot_c "Coolant mass flow rate";
         outer parameter Real Ft = 0.875
            "temperature factor (1 = countercurrent)";                              // minimum value is 0.75

        SI.ThermalResistance R_th "K/W";

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  BED GEOMTERY  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

         outer parameter SI.Length L
            "tank length - used if Calculated_Volume=true";
        // outer parameter SI.Length d_tank
          //  "tank diameter - used if Calculated_Volume=true";
        // outer parameter SI.Volume V_input
          //  "Inner volume of cylindrical tank - used if Calculated_Volume=false";

        // SI.Length d = (if Calculated_Volume == false then (4*V_input/(Modelica.Constants.pi*L))^(1/2) else d_tank) "tank diameter";

        outer SI.Length d;
        outer SI.Volume V_in;

        SI.Volume V_tube_portion;
        // final parameter SI.Length dx = 0.5e-2/3; // thickness of the doscretized mh_bed layer
        outer parameter SI.Length L_i;  // OUTER
        outer parameter SI.Length dx;  // infinitesimal thickness of the external surface (so it is mre representative of the boundary where the Dirichlet condition is defined for the transferred heat to the coolant)

         outer parameter Real  eps;  // porosity (Vgas/V_MH = Vgas/(Vgas+Vsolid))
         outer parameter SI.Length ID;  // Inner tube inner-diameter

        SI.Area A_hex;
        SI.Length Local_Radius;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  INITIAL CONDITIONS  **************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        SI.Mass m_s_tube; // absorbing alloy mass in each tube

         outer parameter SI.Temperature Tamb
            "Initial bed and hydrogen temperature";
        SI.Mass m_g0; // initial hydrogen mass stored in the tank = initial hydrogen mass in the pores (so under the assumption that at time =0 the solid is fully discharged)
        SI.Mass m_g0_full;
        SI.Mass m_des0;

        // added parameters in this test component only in order to define the variable p (i.e. pressure ramp)
        outer parameter SI.Pressure  p0 "Pa";  // Questa diventerá una inner property insieme a quella del tank (che deve essere uguale!) e verrá definita da qualche altra parte dove questo modello e Tank saranno messi dentro
        outer parameter SI.Pressure  p_amb;
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  OPERATIVE VARIABLES  *************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        SI.Temp_K T_bed "Bed and hydrogen temperature";
         outer SI.Pressure p "pressure inside the tank";
        SI.Pressure p_eq "equilibrium pressure [Pa]";

        // Coefficients of the energy balance
        CoefficientOfTemperatureDerivative A
            "coefficient of temperature derivative [J/K]";
        SI.Volume B "coefficient of pressure derivative [m3]";
        SI.Heat E "coefficient of gravimetric density derivative [J]";
        SI.HeatFlowRate D "known term [W]";

        // Masses

        SI.Mass m_g_tot_tube; // total hydrogen mass stored in the tank = absorbed + gaseous in the pores [kg]
        SI.Mass m_gas_tube; // hydrogen mass stored in the pores [kg]
        SI.Mass m_abs_tube; // total hydrogen mass absorbed in the MH [kg]
        SI.MassFlowRate mDot_gas_tube;

        Real F_rp "reaction progress [kg/kg]";

        SI.HeatFlowRate q_gen_p; // Heat of compression [W]
        SI.HeatFlowRate q_gen_DeltaH; // Heat genrated during absorption reaction [W]
        SI.HeatFlowRate q_gen_tot; // Total heat generated in the tank = compession + reaction [W]

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  PROPERTIES  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        // Metal hydride properties stored in records [1]:
         outer replaceable parameter MHSS.Properties.MH.BaseProperties MH_properties
            constrainedby MHSS.Properties.MH.BaseProperties                                    annotation (choicesAllMatching=true);

        GravimetricDensity w;
        SI.Density rho_eff;  // effective density of the MH =
        SI.SpecificHeatCapacity C_MH; // effective thermal capacity = (eps*rho_gas*cp_gas + (1-eps)*rhos_s*cp_s)
        SI.ThermalConductivity k_eff; // effective thermal conductivity = (eps*k_g + (1-eps)*k_s)
         outer parameter SI.SpecificHeatCapacity cp_c; // OUTER
          constant SI.MolarEntropy R = 8.3142 "J/(mol K)";
         constant MolecularWeight MW_g = 2
            "Hydrogen molecular weight in [g/mol]";

         // Assignation of values stored in the record into the variables defined in the present model
        final parameter GravimetricDensity w_max = MH_properties.w_max
            "Max gravimetric density for the material [kgH2/kgTOT]";
        final parameter MolarEnthalpy  DeltaH = MH_properties.DeltaH
            "Enthalpy of generation [J/mol]";
            final parameter MolarEnthalpy  DeltaH_des = MH_properties.DeltaH_des
            "Enthalpy of generation [J/mol]";
        final parameter SI.MolarEntropy   DeltaS = MH_properties.DeltaS
            "Entropy of generation [J/(mol*K)]";
            final parameter SI.MolarEntropy   DeltaS_des = MH_properties.DeltaS_des
            "Entropy of generation [J/(mol*K)]";
        constant SI.Frequency  Ca = MH_properties.Ca "Activation rate [1/s]";
        constant SI.Frequency  Ca_des = MH_properties.Ca_des
            "Activation rate [1/s]";
         constant MolarEnthalpy  Ea = MH_properties.Ea
            "Activation Energy [J/molH2]";
          constant MolarEnthalpy  Ea_des = MH_properties.Ea_des
            "Activation Energy [J/molH2]";
        final parameter SI.Density   rho_s = MH_properties.rho_s
            "Crystalline density [kg/m3]";
        final parameter SI.SpecificHeatCapacityAtConstantPressure  cp_s0 = MH_properties.cp_s0
            "Specific heat capacity [J/(kg K)]";
        final parameter SI.ThermalConductivity  k_s = MH_properties.k_s
            "Thermal conductivity [W/(m K)]";
        /*

final outer parameter GravimetricDensity w_max 
    "Max gravimetric density for the material [kgH2/kgTOT]";
final outer parameter MolarEnthalpy  DeltaH "Enthalpy of generation [J/mol]";
final outer parameter SI.MolarEntropy   DeltaS 
    "Entropy of generation [J/(mol*K)]";
outer parameter SI.Frequency  Ca "Activation rate [1/s]";
 outer parameter MolarEnthalpy  Ea "Activation Energy [J/molH2]";
final outer parameter SI.Density   rho_s "Crystalline density [kg/m3]";
    final outer parameter SI.SpecificHeatCapacityAtConstantPressure  cp_s0 
    "Specific heat capacity [J/(kg K)]";
final outer parameter SI.ThermalConductivity  k_s 
    "Thermal conductivity [W/(m K)]";
*/

        // hydrogen properties (CoolProp):
        replaceable package Medium =
            ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
            Modelica.Media.Interfaces.PartialMedium   annotation (choicesAllMatching=true);

         // Coolprop states and properties
         Medium.ThermodynamicState hydrogen0;
         Medium.ThermodynamicState hydrogen;
         Medium.ThermodynamicState hydrogen1;
         Medium.ThermodynamicState hydrogen2;
         Medium.ThermodynamicState hydrogen3;
         DensityDerivative drhodt;
         Real mu_g;
         SI.SpecificEnthalpy h_g;
         SI.SpecificVolume v_g;
         Real rho_dP;
         SI.SpecificHeatCapacity cp_dP;
         Real rho_dT;
         Real cv_dT;
         Real drhodP_T;
         Real drhodT_P;
         Real dcvdT_P;
         Real dcpdP_T;
         SI.Density rho_g;
         SI.Density rho_g0;
         SI.Density rho_g0_full;
         SI.SpecificInternalEnergy u_g;
         SI.SpecificHeatCapacity  cp_g;
         SI.SpecificEnthalpy h_g_in;
         SI.SpecificEntropy s_g;
         SI.SpecificHeatCapacity  cv_g;
         IsobaricThermalExpansion a_pg;
         IsothermalCompressibility k_tg;
        final parameter Real dP = 1000; // small pressure increment used in the density derivative "drhodt"
        final parameter Real dT = 0.1; //  small temperature increment used in the density derivative "drhodt"

        /*****************************************************************  LINEAR INTERPOLATION FUNCTIONS  ***************************************************************************************************/

           function DATA_interpolation
            input Real x;
            output Real y;
          protected
            Real iNew;
            Real F_data[:]={0,0.02,0.04,0.06,0.13,0.33,0.5,0.66,0.83,0.9};
            Real C_MH_data[:]={500,520,550,555,530,580,560,670,980,1050};
           algorithm
            (
           y,iNew) := Modelica.Math.Vectors.interpolate(
                 F_data,
                 C_MH_data,
                 x);
           end DATA_interpolation;

        Real x_in;
        Real y_out;

         function DATA_interpolation1
            input Real x_k_eff;
            output Real y_k_eff;
          protected
            Real iNew;
            Real p_data[:]={101300, 290000,640000,1480000,3420000,5190000,7060000,10400000,14000000,16400000,17300000,20400000,22900000,25300000};
            Real k_eff_data[:]={0.3070, 0.3077, 0.3848, 0.4606, 0.5281, 0.5572, 0.5769, 0.6030, 0.6207, 0.6534, 0.6973, 0.6672, 0.6737, 0.6677};
         algorithm
            (
           y_k_eff,iNew) := Modelica.Math.Vectors.interpolate(
                 p_data,
                 k_eff_data,
                 x_k_eff);
         end DATA_interpolation1;

        Real x_k_eff_in;
        Real y_k_eff_out;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  INITIAL EQUATIONS  ****************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

           Ports.HeatTransferCoefficientPort_1D B_Port
            annotation (Placement(transformation(extent={{-10,10},{10,30}}),
                iconTransformation(extent={{-10,18},{10,38}})));

          Ports.HeatTransferCoefficientPort_1D A_Port
            annotation (Placement(transformation(extent={{-10,-30},{10,-10}}),
                iconTransformation(extent={{-10,-34},{10,-14}})));
          Ports.Nt_Port Nt_Port
            annotation (Placement(transformation(extent={{-92,-8},{-72,12}}),
                iconTransformation(extent={{-92,-8},{-72,12}})));
          Ports.Mass_Port Mass_Port
            annotation (Placement(transformation(extent={{82,-8},{102,12}})));
        initial equation
          if Charging == true then
        m_g_tot_tube = m_g0;
        w = 0;
        m_abs_tube = 0;
         else
        m_g_tot_tube = m_g0_full + m_des0;
        w = 9/10*w_max;
        m_abs_tube = 9/10*m_des0;
          end if;

        T_bed = Tamb;
          /*
m_g_tot_tube = m_g0;
T_bed = Tamb;
w = 0;
m_gas_tube = m_g0;
m_abs_tube =0;
*/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  EQUATIONS  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        equation

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  EQUILIBRIUM PRESSURE  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        if Charging == true then
          p_eq = p_amb * exp(DeltaH / (R * T_bed) - DeltaS / R); // when T_bed >> => p_eq >> ;  when T_bed <<  =>  p_eq <<
        else
          p_eq = p_amb * exp(DeltaH_des / (R * T_bed) - DeltaS_des / R); // when T_bed >> => p_eq >> ;  when T_bed <<  =>  p_eq <<

        end if;
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  KINETICS EQUATION  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        if Charging == true then
          if p>=p_eq then

             der(w)= Ca * exp(-Ea / (R * T_bed)) * log(p / p_eq) * (w_max - w);

           else

             der(w) = 0;

           end if;
        else
          if p>=p_eq then
            der(F_rp) = 0;
        //      der(w) = 0;

          else
               der(F_rp)= Ca_des * exp(-Ea_des / (R * T_bed)) * log(p / p_eq) * (F_rp);
          //     der(w)= Ca_des * exp(-Ea_des / (R * T_bed)) * log(p / p_eq) * (w_max - w);
          end if;
        end if;

        F_rp = w/w_max; // Reaction progress [2]

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  HYDROGEN MASS FLOW RATE  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        m_s_tube = rho_s*(1-eps)*pi*(ID/2)^2*L;

        m_g0 = rho_g0*V_tube_portion*eps;
        m_g0_full = rho_g0_full*V_tube_portion*eps;
        m_des0 = (1 - eps) * V_tube_portion * rho_s * w_max;

        der(m_g_tot_tube) = eps * V_tube_portion * drhodt +  (1 - eps) * V_tube_portion * rho_s * der(w); // fuelled hydrogen mass flow rate (hydrogen in the pores and hydrogen  in absorbed)

        // der(m_gas_tube) = eps * V_tube_portion * drhodt;
        der(m_abs_tube) = (1 - eps) * V_tube_portion * rho_s * der(w);

        m_gas_tube = eps * V_tube_portion * rho_g;
        mDot_gas_tube = eps * V_tube_portion * drhodt;

        Mass_Port.m_tot = m_g_tot_tube;

        /*
m_g0 = rho_g0*V_tube_portion*eps;
der(m_g_tot_tube) = eps * V_tube_portion * drhodt +  (1 - eps) * V_tube_portion * rho_s * der(w); // fuelled hydrogen mass flow rate (hydrogen in the pores and hydrogen in absorbed)

der(m_gas_tube) = eps * V_tube_portion * drhodt;
der(m_abs_tube) = (1 - eps) * V_tube_portion * rho_s * der(w);

Mass_Port.m_tot = m_g_tot_tube;
// Mass_Port.m_g_tot = m_gas_tube;
// Mass_Port.m_abs_tot = m_abs_tube;
*/
        /***************************************************************************************************************************************************************
******************************************************************  SOLID BED TEMPERATURE  *********************************************************************
****************************************************************************************************************************************************************/

         A = ((1-eps)*rho_s*w*cp_g + rho_eff*C_MH)*V_tube_portion;  // Coefficient of Temperature derivative - Considers the effective heat capacity C_MH = eps*rho_g*cp_g + (1-eps)*rho_s*cp_s = 500 J/(kg m3)

         B = ((1-eps) + ((1-eps)*rho_s/rho_g*w+eps)*(1 - a_pg*T_bed) - eps)*V_tube_portion; // Coefficient of pressure derivative

        if Charging == true then
        D = (der(m_g_tot_tube)*(h_g_in - h_g) + A_Port.qDot + B_Port.qDot); // known term (energy in, and its time variation in the tank, and exchanged heat)
         // D = (der(m_g_tot_tube)*(h_g_in-h_g) + A_Port.qDot + B_Port.qDot); // known term (energy in, and its time variation in the tank, and exchanged heat)

         E = ((1-eps)*rho_s * DeltaH/(MW_g/1000))*V_tube_portion;  // Coefficient of gravimetric-density derivative

        else
          D = (der(m_g_tot_tube)*(h_g_in - h_g) + A_Port.qDot + B_Port.qDot);    // known term (energy in, and its time variation in the tank, and exchanged heat)
          E = ((1-eps)*rho_s * DeltaH_des/(MW_g/1000))*V_tube_portion;  // Coefficient of gravimetric-density derivative

        end if;

        // E = ((1-eps)*rho_s * DeltaH/(MW_g/1000))*V_tube_portion;  // Coefficient of gravimetric-density derivative

         // D = (der(m_g_tot_tube)*(h_g_in-h_g) + A_Port.qDot + B_Port.qDot); // known term (energy in, and its time variation in the tank, and exchanged heat)

        A_hex = 2*pi*(ID/2-(A_Port.Counter-1)*dx)*L_i*Nt_Port.Nt;
        V_tube_portion = pi*((ID/2-(A_Port.Counter-1)*dx)^2 - (ID/2-A_Port.Counter*dx)^2)*L_i*Nt_Port.Nt;

        Local_Radius = ID/2 - A_Port.Counter*dx;
        R_th = dx/(A_hex*k_eff); // thermal resistance

        A_Port.qDot = (A_Port.T-T_bed)/(R_th);
        B_Port.qDot = (B_Port.T-T_bed)/(R_th); // the effect of the number of tubes is already included in the definition of A_hex

        B_Port.Counter = A_Port.Counter + 1;

        der(T_bed)= (D - B*der(p) - E *der(w))/A; // Energy balance

        /***************************************************************************************************************************************************************
  ***************************************************************  MH properties  ******************************************************************************
  **************************************************************************************************************************************************************/
        rho_eff = (1-eps)*rho_s + eps * rho_g; // effective density

        x_in = F_rp;
        y_out = DATA_interpolation(x=x_in);

        x_k_eff_in = p;
        y_k_eff_out = DATA_interpolation1(x_k_eff=x_k_eff_in);

        if VariableProperties == true then

        C_MH = y_out;
        k_eff = y_k_eff_out;

        else
        k_eff = MH_properties.k_s;
          C_MH = 500;
        end if;

        /***************************************************************************************************************************************************************
  ***************************************************************  Heat in the tank  ******************************************************************************
  **************************************************************************************************************************************************************/

        q_gen_p = eps*V_tube_portion*der(p); // Heat of compression [1]
        q_gen_DeltaH = -E*der(w); // Heat of reaction
        q_gen_tot = q_gen_p + q_gen_DeltaH;

        /***************************************************************************************************************************************************************
  *******************************************  Thermodynamic states and properties calculations (COOLPROP)  ****************************************************
  **************************************************************************************************************************************************************/

         drhodt = drhodP_T * der(p) + drhodT_P * der(T_bed); // derivative of the hydrogen density

          v_g=1/rho_g;

         hydrogen0=Medium.setState_pT(p=p0, T=Tamb);
         hydrogen=Medium.setState_pT(p=p, T=T_bed);
         hydrogen1=Medium.setState_pT(p=p, T=Tamb);
         hydrogen2=Medium.setState_pT(p=p+dP, T=T_bed);
         hydrogen3=Medium.setState_pT(p=p, T=T_bed+dT);

         h_g_in=Medium.specificEnthalpy(hydrogen1);
         h_g=Medium.specificEnthalpy(hydrogen) "Dynamic enthalpy";
         cp_g=Medium.specificHeatCapacityCp(hydrogen);
         u_g=Medium.specificInternalEnergy(hydrogen);
          rho_g0=Medium.density(hydrogen0);
         rho_g=Medium.density(hydrogen);
           rho_g0_full=Medium.density(hydrogen0);
         s_g=Medium.specificEntropy(hydrogen);
         cv_g=Medium.specificHeatCapacityCv(hydrogen);
         mu_g=Medium.dynamicViscosity(hydrogen);

         rho_dP =Medium.density(hydrogen2);
         cp_dP= Medium.specificHeatCapacityCp(hydrogen2);

         drhodP_T = (rho_dP - rho_g) / dP;
         drhodT_P = (rho_dT - rho_g) / dT;

         a_pg = rho_g*(1/rho_dT - 1/rho_g) / dT;
         k_tg = - rho_g*(1/rho_dP - 1/rho_g)/dP;

         rho_dT = Medium.density(hydrogen3);
         cv_dT = Medium.specificHeatCapacityCv(hydrogen3);

         dcpdP_T = (cp_dP - cp_g) / dP;
         dcvdT_P = (cv_dT - cv_g) / dT;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************ REFERENCES ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /*
[1] M. Visaria, I. Mudawar, T. Pourpoint, and S. Kumar, “Study of heat transfer and kinetics parameters influencing 
    the design of heat exchangers for hydrogen storage in high-pressure metal hydrides,” Int. J. Heat Mass Transf., vol. 
    53, no. 9–10, pp. 2229–2239, Apr. 2010.
[2] T. L. Pourpoint, V. Velagapudi, I. Mudawar, Y. Zheng, and T. S. Fisher, “Active cooling of a metal hydride 
    system for hydrogen storage,” Int. J. Heat Mass Transf., vol. 53, no. 7–8, pp. 1326–1332, Mar. 2010.
[3]
[4]

*/

         annotation (Diagram(graphics), Icon(graphics={Bitmap(extent={{-70,56},{84,-52}},
                    fileName="modelica://Complete_Refuelling/../../Figures for components/MH_bed.png")}));
        end MH_bed;

        block Stored_Mass
          "Partial block with a BooleanVectorInput and a BooleanOutput signal"

          import SI = Modelica.SIunits;

          // extends Modelica.Blocks.Interfaces.MISO(nin=N);
        //  parameter Real k[N]=ones(N) "Optional: sum coefficients";
          outer parameter Integer N;
          SI.Mass m_stored;

        // parameter Integer nin=1 "Number of inputs";

          Ports.Mass_Port mass_Port[N]
            annotation (Placement(transformation(extent={{-122,-18},{-92,18}})));
        Real v[N];
        equation
          for i in 1:N loop
            v[i] = mass_Port[i].m_tot;
          end for;

          m_stored = sum(v);

         // y = k*u;
          // m_stored = y;
          annotation (defaultComponentName="sum1",
            Documentation(info="
<HTML>
<p>
This blocks computes output <b>y</b> as
<i>sum</i> of the elements of the input signal vector
<b>u</b>:
</p>
<pre>
    <b>y</b> = <b>u</b>[1] + <b>u</b>[2] + ...;
</pre>
<p>
Example:
</p>
<pre>
     parameter:   nin = 3;

  results in the following equations:

     y = u[1] + u[2] + u[3];
</pre>

</HTML>
"),         Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            grid={2,2}), graphics={
                               Text(
              extent={{-150,150},{150,110}},
              textString="%name",
              lineColor={0,0,255}),
                Rectangle(
                  extent={{-66,62},{60,-68}},
                  lineColor={0,0,0},
                  fillColor={170,85,255},
                  fillPattern=FillPattern.Solid),
                                   Line(
              points={{24,36},{-36,36},{4,-4},{-36,-44},{24,-44}},
              color={0,0,0},
              thickness=1)}),
            Diagram(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            grid={2,2}), graphics={Rectangle(
              extent={{-100,-100},{100,100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid), Line(
              points={{26,42},{-34,42},{6,2},{-34,-38},{26,-38}},
              color={0,0,0},
              thickness=0.25)}));
        end Stored_Mass;

        model Nt_distributor

          Ports.Nt_Port Nt_Port_out
            annotation (Placement(transformation(extent={{-12,72},{8,92}})));
        //Ports.Nt_Port Nt_Port_out2
        //    annotation (Placement(transformation(extent={{44,-10},{64,10}})));

          Ports.Nt_Port Nt_Port_in
            annotation (Placement(transformation(extent={{-10,-106},{10,-86}})));

        equation
        Nt_Port_in.Nt = Nt_Port_out.Nt;
        //Nt_Port_in.Nt = Nt_Port_out2.Nt;

        annotation (Icon(graphics={Rectangle(
                  extent={{-56,78},{66,-92}},
                  lineColor={0,0,0},
                  fillColor={255,128,0},
                  fillPattern=FillPattern.Solid), Text(
                  extent={{-34,38},{40,-50}},
                  lineColor={0,0,0},
                  lineThickness=1,
                  fillColor={170,85,255},
                  fillPattern=FillPattern.Solid,
                  textString="Nt")}),               Diagram(graphics));
        end Nt_distributor;

        model External_Coolant

          import SI = Modelica.SIunits;

          Ports.Temperature_port C_in annotation (Placement(transformation(extent={{-80,
                    -4},{-60,16}}), iconTransformation(extent={{-92,-14},{-60,16}})));
          Ports.Temperature_port C_out annotation (Placement(transformation(extent={{74,
                    -6},{94,14}}), iconTransformation(extent={{62,-16},{94,14}})));

        outer parameter SI.MassFlowRate mDot_c "Coolant mass flow rate";
        outer parameter SI.SpecificHeatCapacity cp_c;

          Ports.HeatTransferCoefficientPort_1D heatTransferCoefficientPort_1D
            annotation (Placement(transformation(extent={{-10,54},{10,74}})));
        equation
          // q_hex = heatTransferCoefficientPort_EXT_1D.q must be the total heat transferred (i.e. related to ALL the tubes)

        heatTransferCoefficientPort_1D.T = (C_out.T  +  C_in.T)/2;
        C_out.T = heatTransferCoefficientPort_1D.qDot/(mDot_c*cp_c)  +  C_in.T;

          annotation (Icon(graphics={Rectangle(
                  extent={{-72,60},{78,-62}},
                  lineColor={0,0,255},
                  fillColor={85,170,255},
                  fillPattern=FillPattern.Solid)}), Diagram(graphics));
        end External_Coolant;
      end One_Dimensional_Models;

      package Zero_Dimensional_Models
        model Tank_system_with_Hexs
         // import HPMH_Refuelling_Paper1;

          import SI = Modelica.SIunits;

        type GravimetricDensity = Real(quantity = "Gravimetric Density", unit = "kg/kg");
        type MolarEnthalpy = Real(quantity = "MolarEnthalpy", unit = "J/mol");
        type GravimetricDensity_syst = Real(quantity = "Gravimetric Density", unit = "kg_H2/kg_syst");
        type VolumetricDensity_syst = Real(quantity = "Gravimetric Density", unit = "kg_H2/kg_syst");

        /*******************************  PORTS  **************************************************************************************************/

         Ports.p_T p_T
            annotation (Placement(transformation(extent={{-10,10},{10,30}}),
                iconTransformation(extent={{80,-10},{100,10}})),Diagram(graphics));
         Ports.PressurePort p_ramp
            annotation (Placement(transformation(extent={{-88,-10},{-68,10}}),
                iconTransformation(extent={{-88,-10},{-68,10}})));
        /******************************  GRAPHICS  ***************************************************************************************************/

        /******************************  BOOLEAN  ***************************************************************************************************/
         inner parameter Boolean Calculated_Volume = true
            "if true then the volume is calculated from L_input and d_tank";
         inner parameter Boolean VariableProperties = true
            "if true then C_MH and k_MH are interpolated from experimental values (p,Reation_Progress)";
         inner parameter Boolean Bell_Delaware_Method = true
            "If false then Kern Method";
        inner parameter Boolean Charging = true "if false then discharging";

            /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  Volumetric and Gravimetric densities definition  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        GravimetricDensity_syst Grav_Density;
        VolumetricDensity_syst Vol_Density; // it is assumed that the outer volume of the tank coincides with the inenr volume, as for a tubular tank design the vessel does not have to bear the high pressure (in the tubes) and can be made of a "thin" plastic
        SI.Mass HEX_weight;
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  COOLANT AND HEX ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        inner parameter SI.MassFlowRate mDot_c = 30 "Coolant mass flow rate";
        inner parameter Real Ft = 0.875
            "temperature factor (1 = countercurrent)";                             // minimum value is 0.75

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  BED GEOMTERY  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        inner parameter Real eps=0.6 "Bed porosity [m3_g/m3_MH]"; // outer in MH_bed

         inner parameter SI.Length L = 1
            "tank length - used if Calculated_Volume=true";
         final inner parameter SI.Length Ltube = L; // Then assumption that the tube length (HEX) equals the tank length
         inner parameter SI.Length d_tank = 0.4;
        // final inner parameter SI.Length Ds = d_tank
         //   "tank diameter - used if Calculated_Volume=true";                                            // The tank inner parameter equals the inner shell diameter of the Shell and Tube
        inner parameter SI.Volume V_input = 0.180
            "Inner volume of cylindrical tank - used if Calculated_Volume=false";
        final inner parameter SI.Length ID = 2*delta
            "inner tube diameter = 2*delta";
        final inner parameter SI.Length OD = ID/0.8;
         parameter SI.Length delta = 0.5e-2
            "Critical MH thickness = 1/2 inner tube diameter";

        inner SI.Length d = (if Calculated_Volume == false then (4*V_input/(Modelica.Constants.pi*L))^(1/2) else d_tank)
            "tank diameter";

        inner SI.Volume V_in = (if Calculated_Volume == false then V_input else Modelica.Constants.pi*(d^2/4)*L)
            "Inner volume of the cylindrical tank";

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  INITIAL CONDITIONS  **************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        inner SI.Pressure p0 "Initial tank pressure";
        inner SI.Temperature Tamb_start "Bed temperature";
        inner SI.Pressure p_amb "Ambient pressure";

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  PRESSURE RAMP  **************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        inner SI.Pressure p; // =p_ramp

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  MH AND COOLANT PROPERTIES  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        // Metal hydride properties stored in records [1]:
        inner replaceable parameter MHSS.Properties.MH.BaseProperties MH_properties
            constrainedby MHSS.Properties.MH.BaseProperties                     annotation (choicesAllMatching=true);

        // Assignation of values stored in the record into the variables defined in the present model
        final inner parameter GravimetricDensity w_max =  MH_properties.w_max
            "Max gravimetric density for the material [kgH2/kgTOT]";
        final inner parameter MolarEnthalpy DeltaH = MH_properties.DeltaH
            "Enthalpy of generation [J/mol]";
            final inner parameter MolarEnthalpy  DeltaH_des = MH_properties.DeltaH_des
            "Enthalpy of generation [J/mol]";
        final inner parameter SI.MolarEntropy DeltaS = MH_properties.DeltaS
            "Entropy of generation [J/(mol*K)]";
                final inner parameter SI.MolarEntropy   DeltaS_des = MH_properties.DeltaS_des
            "Entropy of generation [J/(mol*K)]";
        inner constant SI.Frequency Ca = MH_properties.Ca
            "Activation rate [1/s]";
        inner constant SI.Frequency  Ca_des = MH_properties.Ca_des
            "Activation rate [1/s]";
        inner constant MolarEnthalpy Ea = MH_properties.Ea
            "Activation Energy [J/molH2]";
            inner constant MolarEnthalpy  Ea_des = MH_properties.Ea_des
            "Activation Energy [J/molH2]";

        final inner parameter SI.Density rho_s = MH_properties.rho_s
            "Crystalline density [kg/m3]";
        final inner parameter SI.SpecificHeatCapacityAtConstantPressure cp_s0= MH_properties.cp_s0
            "Specific heat capacity [J/(kg K)]";
        final inner parameter SI.ThermalConductivity k_s = MH_properties.k_s
            "Thermal conductivity [W/(m K)]";

        inner replaceable parameter
            MHSS.Properties.Coolant.CoolantBaseProperties Coolant
            constrainedby MHSS.Properties.Coolant.CoolantBaseProperties                                         annotation (choicesAllMatching=true);

        final inner parameter SI.Density rho_c = Coolant.rho_c "kg/m3";
        final inner parameter SI.SpecificHeatCapacity cp_c = Coolant.cp_c
            "J/(kg K)";
        final inner parameter SI.ThermalConductivity k_c = Coolant.k_c "W/(m K)";
        final inner parameter SI.DynamicViscosity mu_c = Coolant.mu_c "Pa s";

            // Tube properties stored in records:
        inner replaceable parameter MHSS.Properties.Tube.BaseProperties Tube_properties
            constrainedby MHSS.Properties.Tube.BaseProperties                                    annotation (choicesAllMatching=true);

            final parameter SI.ThermalConductivity  k_tube = Tube_properties.k_tube;
            final parameter SI.Density  rho_tube = Tube_properties.rho_tube;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  HEAT EXCHANGER  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        inner replaceable parameter
            MHSS.Heat_Transfer.Heat_Exchangers.Heat_Exchanger_Configuration HEXCHANGER_USED
            constrainedby
            MHSS.Heat_Transfer.Heat_Exchangers.Heat_Exchanger_Configuration                annotation (choicesAllMatching=true);

         MHSS.Heat_Transfer.Heat_Exchangers.Heat_exchangers heat_exchangers(
              HEXCHANGER=HEXCHANGER_USED)                                                                              annotation (Placement(transformation(extent={{-52,
                    -120},{36,-8}})));

        inner parameter SI.Temperature Tc_in = 273.15
            "Inlet coolant temperature";
        inner parameter Real PR = 1.25
            "Tube pitch to tube diameter ratio: Pt/OD";                             // Tube pitch to tube outer diameter ratio (Pt/OD) http://www.engineeringpage.com/technology/thermal/pitch.html
        inner parameter Modelica.SIunits.Conversions.NonSIunits.Angle_deg theta_tp = 30
            "Layout angle (30,60,45,90deg)";                                                                        // in deg
        inner parameter Integer Nt_input = 100
            "Used if NumberOfTubes_Given = true";
        inner parameter Real Nss = 1
            "Number of Sealing strips (pairs) in one baffle spacing";
        inner parameter Real B0=0.5
            "Initial value for the baffle space 0<B<1 - percentage (>0.2 (>5cm in general))";                       // Baffle space as percentage of the shell diameter "Ds"
        inner parameter Real Bc = 25e-2
            "Baffle cut (20-49%) with 20-25% optimum";
        inner parameter Integer Nb = 4
            "Number of baffles - input if NumberOfTubes_Given = true";
        inner parameter Real Lo_star = 1
            "1=<Lo_star<2 --> Ratio between the inlet baffle spacing and the central baffle spacing"; // Ratio between the outlet baffle spacing and the central baffle spacing
        inner parameter Real Li_star = 1
            "1=<Li_star<2 --> Ratio between the inlet baffle spacing and the central baffle spacing"; // Ratio between the inlet baffle spacing and the central baffle spacing
        inner parameter Integer N_shell = 1 "Number of shells";
        inner parameter Real a_r = 0.75 "Aspect ratio; 0<a<1 (ID/OD)";

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  EQUATIONS  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        //  Heat_Transfer.Heat_Exchangers.Heat_exchangers heat_exchangers
          MH_bed mH_bed
            annotation (Placement(transformation(extent={{-92,-66},{76,60}})));
        equation
         p_ramp.p = p; // p_ramp coming from Pressure_ramp_and_Tamb

         p_T.T = Tamb_start; // ambient temperature
         p_T.p0 = p0;   // initial pressure in the tank
         p_T.p_amb = p_amb;

        Grav_Density = mH_bed.m_g_tot/(mH_bed.m_g_tot + mH_bed.m_s_tot + HEX_weight); // it is assumed that the outer volume of the tank coincides with the inenr volume, as for a tubular tank design the vessel does not have to bear the high pressure (in the tubes) and can be made of a "thin" plastic
        Vol_Density = mH_bed.m_g_tot/(V_in*1e3);
        HEX_weight = rho_tube*mH_bed.heatTransferCoefficientPort.Nt*Modelica.Constants.pi*((OD/2)^2-(ID/2)^2)*L;

          connect(mH_bed.heatTransferCoefficientPort, heat_exchangers.heatTransferCoefficientPort)
            annotation (Line(
              points={{-8,-18.12},{-8,-41.6}},
              color={0,0,0},
              smooth=Smooth.None));
            annotation (Placement(transformation(extent={{-30,-98},{40,-26}})),
                     Icon(graphics={                                                                                                    Text(
                  extent={{-20,60},{30,26}},
                  lineColor={0,0,0},
                  fillColor={170,85,255},
                  fillPattern=FillPattern.Solid,
                  textString="0D"), Bitmap(extent={{-68,-50},{84,54}}, fileName=
                     "modelica://HySDeP/Graphics/Tank1.jpg")}),
              Diagram(graphics),
                        Placement(transformation(extent={{-94,-52},{100,42}})));
        end Tank_system_with_Hexs;

        model MH_bed

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  UNITS    ***************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

          import SI = Modelica.SIunits;
          import Modelica.Constants.pi;

        type GravimetricDensity = Real(quantity = "Gravimetric Density", unit = "kg/kg");
        type MolarEnthalpy = Real(quantity = "MolarEnthalpy", unit = "J/mol");
        type MolecularWeight = Real(quantity= "MolecularWeight", unit="g/mol");
        type CoefficientOfTemperatureDerivative = Real(quantity= "CoefficientOfTemperatureDerivative", unit="J/K");
        type IsobaricThermalExpansion = Real(quantity="IsobaricThermalExpansion", unit="1/K");
        type IsothermalCompressibility = Real(quantity="IsobaricThermalExpansion", unit="1/Pa");
        type DensityDerivative = Real(quantity="DensityDerivative", unit="(kg/m3)/s");
        type DensityDerivativeAtConstantTemperature = Real(quantity="DensityDerivativeAtConstantTemperature", unit="(kg/m3)/Pa");
        type DensityDerivativeAtConstantPressure = Real(quantity="DensityDerivativeAtConstantPressure", unit="(kg/m3)/K");
        type ThermalContactResistance = Real(quantity="ThermalContactResistance", unit="(m2.K)/W");

        /************************************************  GRAPHICS  ************************************************************************************************************************************************************************************************************************************/

          Ports.HeatTransferCoefficientPort heatTransferCoefficientPort annotation (Placement(transformation(extent={{-14,-36},
                    {12,-12}}),                                                                                                    iconTransformation(extent={{-12,-36},{12,-12}})));

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  BOOLEAN  *************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        outer parameter Boolean Calculated_Volume
            "if true then the volume is calculated from L_input and d_tank";

        outer parameter Boolean VariableProperties
            "if true then C_MH and k_MH are interpolated from experimental values (p,Reation_Progress)";
        outer parameter Boolean Charging "if false then discharging";

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  COOLANT AND HEX **************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        outer parameter SI.MassFlowRate mDot_c = 30 "Coolant mass flow rate";
        outer parameter Real Ft = 0.875
            "temperature factor (1 = countercurrent)";                             // minimum value is 0.75
         SI.HeatFlowRate q_hex;

        SI.CoefficientOfHeatTransfer h;
        SI.CoefficientOfHeatTransfer U;
        parameter ThermalContactResistance R_th_cont = 2000e-6 "[(m2.K)/W]"; //  contact resistance between tube/coolant and hydride [1]
        ThermalContactResistance R_th_tube "[(m2.K)/W]";
        /*
SI.TemperatureDifference DeltaT_lm;
SI.TemperatureDifference DeltaT_lm_NUM(start=20);
SI.TemperatureDifference DeltaT_lm_DEN;
// SI.TemperatureDifference DeltaT_lm_DEN_0;

SI.TemperatureDifference DeltaT_lm_DEN1;
SI.TemperatureDifference DeltaT_lm_DEN2;
*/
        outer parameter SI.Temperature Tc_in;
        outer parameter SI.Pressure p_amb;

        // final parameter SI.Temperature Tc_out = 273.15+1e-5;
        SI.Temperature Tc_out(start=Tc_in); // the term "1e-3" is added to avoid division by 0 in DeltaT_lm at around time 0
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  BED GEOMTERY  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

         outer parameter SI.Length L
            "tank length - used if Calculated_Volume=true";
        // outer parameter SI.Length d_tank     "tank diameter - used if Calculated_Volume=true";
            // outer parameter SI.Volume V_input   "Inner volume of cylindrical tank - used if Calculated_Volume=false";

        outer parameter SI.Length d;
        //= (if Calculated_Volume == false then (4*V_input/(Modelica.Constants.pi*L))^(1/2) else d_tank) "tank diameter";

        outer parameter SI.Volume V_in;
        //= (if Calculated_Volume == false then V_input else Modelica.Constants.pi*(d^2/4)*L)     "Inner volume of the cylindrical tank";

        SI.Volume V_MH_tot;
        final parameter SI.Volume V_tube = pi*ID^2/4*L;

        outer parameter Real  eps; // porosity (Vgas/V_MH = Vgas/(Vgas+Vsolid))
        outer parameter SI.Length ID = 2*5e-3; // Inner tube inner-diameter
        outer parameter SI.Length OD = ID/0.8; // Inner tube outer-diameter

        final parameter SI.Area A_in = pi*ID*L;
        final parameter SI.Area A_ext = pi*OD*L;
        Real V_ratio; // Ratio between MH volume (in tubes) and inner tank volume
        Real V_ratio_OD; // Ratio between MH volume (in tubes) and inner tank volume

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  INITIAL CONDITIONS  **************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        outer parameter SI.Pressure  p0 "Initial hydrogen pressure in the tank";  // Questa diventerá una outer property insieme a quella definita nella source (che deve essere uguale!) e verrá definita da qualche altra parte dove questo modello e Tank saranno messi dentro
        outer parameter SI.Temperature  Tamb_start
            "Initial bed and hydrogen temperature";
        SI.Mass m_g0; // initial hydrogen mass stored in the tank = initial hydrogen mass in the pores (so under the assumption that at time =0 the solid is fully discharged)
        SI.Mass m_g0_full;
        SI.Mass m_des0;
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  OPERATIVE VARIABLES  *************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        SI.Temp_K T_bed "Bed and hydrogen temperature";
        outer SI.Pressure p "pressure inside the tank";
        SI.Pressure p_eq "equilibrium pressure [Pa]";

        // Coefficients of the energy balance
        CoefficientOfTemperatureDerivative A
            "coefficient of temperature derivative [J/K]";
        SI.Volume B "coefficient of pressure derivative [m3]";
        SI.Heat E "coefficient of gravimetric density derivative [J]";
        SI.HeatFlowRate D "known term [W]";

        // Masses
        SI.Mass m_g_tot(min=0); // total hydrogen mass stored in the tank = absorbed + gaseous in the pores [kg]
        SI.Mass m_gas(min=0); // hydrogen mass stored in the pores [kg]
        SI.Mass m_abs(min=0); // total hydrogen mass absorbed in the MH [kg]
        SI.Mass m_g_tot_tube(min=0); // total hydrogen mass stored in the tank = absorbed + gaseous in the pores [kg]
        SI.Mass m_gas_tube(min=0); // hydrogen mass stored in the pores [kg]
        SI.Mass m_abs_tube(min=0); // total hydrogen mass absorbed in the MH [kg]
        SI.Mass m_s_tot(min=0); // mass of the absorbent in the tank (i.e. MH) [kg]
        SI.Mass m_s_tube(min=0); // mass of the absorbent per tube (i.e. MH) [kg]
        SI.MassFlowRate mDot_gas_tube;

        Real F_rp(min=0, max=1) "reaction progress [kg/kg]";
        SI.HeatFlowRate q_gen_p_tube; // Heat of compression [W]
        SI.HeatFlowRate q_gen_DeltaH_tube; // Heat genrated during absorption reaction [W]
        SI.HeatFlowRate q_gen_tot_tube; // Total heat generated in the tank = compession + reaction [W]
        SI.HeatFlowRate q_gen_p_tot; // Heat of compression [W]
        SI.HeatFlowRate q_gen_DeltaH_tot; // Heat genrated during absorption reaction [W]
        SI.HeatFlowRate q_gen_tot; // Total heat generated in the tank = compession + reaction [W]
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  PROPERTIES  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        // Metal hydride properties stored in records [1]:
        outer replaceable parameter MHSS.Properties.MH.BaseProperties MH_properties
            constrainedby MHSS.Properties.MH.BaseProperties                                    annotation (choicesAllMatching=true);

        GravimetricDensity w(min=0, max=w_max);
        SI.Density rho_eff;  // effective density of the MH =
        SI.SpecificHeatCapacity C_MH; // effective thermal capacity = (eps*rho_gas*cp_gas + (1-eps)*rhos_s*cp_s)
        SI.ThermalConductivity k_eff; // effective thermal conductivity = (eps*k_g + (1-eps)*k_s)
        outer parameter SI.SpecificHeatCapacity cp_c;
         // Assignation of values stored in the record into the variables defined in the present model
        final parameter GravimetricDensity w_max = MH_properties.w_max
            "Max gravimetric density for the material [kgH2/kgTOT]";
        final parameter MolarEnthalpy  DeltaH = MH_properties.DeltaH
            "Enthalpy of generation [J/mol]";
            final parameter MolarEnthalpy  DeltaH_des = MH_properties.DeltaH_des
            "Enthalpy of generation [J/mol]";
        final parameter SI.MolarEntropy   DeltaS = MH_properties.DeltaS
            "Entropy of generation [J/(mol*K)]";
            final parameter SI.MolarEntropy   DeltaS_des = MH_properties.DeltaS_des
            "Entropy of generation [J/(mol*K)]";
        constant SI.Frequency  Ca = MH_properties.Ca "Activation rate [1/s]";
        constant SI.Frequency  Ca_des = MH_properties.Ca_des
            "Activation rate [1/s]";
         constant MolarEnthalpy  Ea = MH_properties.Ea
            "Activation Energy [J/molH2]";
          constant MolarEnthalpy  Ea_des = MH_properties.Ea_des
            "Activation Energy [J/molH2]";
        final parameter SI.Density   rho_s = MH_properties.rho_s
            "Crystalline density [kg/m3]";
        final parameter SI.SpecificHeatCapacityAtConstantPressure  cp_s0 = MH_properties.cp_s0
            "Specific heat capacity [J/(kg K)]";
        final parameter SI.ThermalConductivity  k_s= MH_properties.k_s
            "Thermal conductivity [W/(m K)]";
         constant SI.MolarEntropy R = 8.3142 "J/(mol K)";
         constant MolecularWeight MW_g = 2
            "Hydrogen molecular weight in [g/mol]";

        // hydrogen properties (CoolProp):
        replaceable package Medium =
            ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
            Modelica.Media.Interfaces.PartialMedium   annotation (choicesAllMatching=true);

            // Tube properties stored in records:
        outer replaceable parameter MHSS.Properties.Tube.BaseProperties Tube_properties
            constrainedby MHSS.Properties.Tube.BaseProperties                                    annotation (choicesAllMatching=true);

            final parameter SI.ThermalConductivity  k_tube = Tube_properties.k_tube;

         // Coolprop states and properties
         Medium.ThermodynamicState hydrogen0;
         Medium.ThermodynamicState hydrogen;
         Medium.ThermodynamicState hydrogen1;
         Medium.ThermodynamicState hydrogen2;
         Medium.ThermodynamicState hydrogen3;
         Medium.ThermodynamicState hydrogen4;
         DensityDerivative drhodt;
         Real mu_g;
         SI.SpecificEnthalpy h_g;
         SI.SpecificVolume v_g;
         Real rho_dP;
         SI.SpecificHeatCapacity cp_dP;
         Real rho_dT;
         Real cv_dT;
         Real drhodP_T;
         Real drhodT_P;
         Real dcvdT_P;
         Real dcpdP_T;
         SI.Density rho_g;
         SI.Density rho_g0;
          SI.Density rho_g0_full;
         SI.SpecificInternalEnergy u_g;
         SI.SpecificHeatCapacity  cp_g;
         SI.SpecificEnthalpy h_g_in;
         SI.SpecificEntropy s_g;
         SI.SpecificHeatCapacity  cv_g;
         IsobaricThermalExpansion a_pg;
         IsothermalCompressibility k_tg;
        final parameter Real dP = 1000; // small pressure increment used in the density derivative "drhodt"
        final parameter Real dT = 0.1; //  small temperature increment used in the density derivative "drhodt"

        /*****************************************************************  LINEAR INTERPOLATION FUNCTIONS  ***************************************************************************************************/

           function DATA_interpolation
            input Real x;
            output Real y;
          protected
            Real iNew;
            Real F_data[:]={0,0.02,0.04,0.06,0.13,0.33,0.5,0.66,0.83,0.9};
            Real C_MH_data[:]={500,520,550,555,530,580,560,670,980,1050};
           algorithm
            (
           y,iNew) := Modelica.Math.Vectors.interpolate(
                 F_data,
                 C_MH_data,
                 x);
           end DATA_interpolation;

        Real x_in;
        Real y_out;

         function DATA_interpolation1
            input Real x_k_eff;
            output Real y_k_eff;
          protected
            Real iNew;
            Real p_data[:]={101300, 290000,640000,1480000,3420000,5190000,7060000,10400000,14000000,16400000,17300000,20400000,22900000,25300000};
            Real k_eff_data[:]={0.3070, 0.3077, 0.3848, 0.4606, 0.5281, 0.5572, 0.5769, 0.6030, 0.6207, 0.6534, 0.6973, 0.6672, 0.6737, 0.6677};
         algorithm
            (
           y_k_eff,iNew) := Modelica.Math.Vectors.interpolate(
                 p_data,
                 k_eff_data,
                 x_k_eff);
         end DATA_interpolation1;

        Real x_k_eff_in;
        Real y_k_eff_out;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  INITIAL EQUATIONS  ****************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        initial equation
         if Charging == true then
        m_g_tot_tube = m_g0;
        w = 0;
        m_abs_tube = 0;
         else
        m_g_tot_tube = m_g0_full + m_des0;
        w = 9/10*w_max;
        m_abs_tube = 9/10*m_des0;
         end if;

        T_bed = Tamb_start;
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  EQUATIONS  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        equation

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  EQUILIBRIUM PRESSURE  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        if Charging == true then
          p_eq = p_amb * exp(DeltaH / (R * T_bed) - DeltaS / R); // when T_bed >> => p_eq >> ;  when T_bed <<  =>  p_eq <<
        else
          p_eq = p_amb * exp(DeltaH_des / (R * T_bed) - DeltaS_des / R); // when T_bed >> => p_eq >> ;  when T_bed <<  =>  p_eq <<

        end if;
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  KINETICS EQUATION  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        if Charging == true then
          if p>=p_eq then

             der(w)= Ca * exp(-Ea / (R * T_bed)) * log(p / p_eq) * (w_max - w);

           else

             der(w) = 0;

           end if;
        else
          if p>=p_eq then
            der(F_rp) = 0;
        //      der(w) = 0;

          else
               der(F_rp)= Ca_des * exp(-Ea_des / (R * T_bed)) * log(p / p_eq) * (F_rp);
          //     der(w)= Ca_des * exp(-Ea_des / (R * T_bed)) * log(p / p_eq) * (w_max - w);
          end if;
        end if;

        F_rp = w/w_max; // Reaction progress [2]

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  HYDROGEN MASS FLOW RATE  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        m_g0 = rho_g0*V_tube*eps;
        m_g0_full = rho_g0_full*V_tube*eps;
        m_des0 = (1 - eps) * V_tube * rho_s * w_max;

        der(m_g_tot_tube) = eps * V_tube * drhodt +  (1 - eps) * V_tube * rho_s * der(w); // fuelled hydrogen mass flow rate (hydrogen in the pores and hydrogen  in absorbed)

        // der(m_gas_tube) = eps * V_tube * drhodt;
        der(m_abs_tube) = (1 - eps) * V_tube * rho_s * der(w);

        m_gas_tube = eps * V_tube * rho_g;
        mDot_gas_tube = eps * V_tube * drhodt;
        m_g_tot = m_g_tot_tube*heatTransferCoefficientPort.Nt;
        m_gas = m_gas_tube*heatTransferCoefficientPort.Nt;
        m_abs = m_abs_tube*heatTransferCoefficientPort.Nt;

        // mass of absorbent (i.e. MH) in the tank
        m_s_tube = (1-eps)*V_tube*rho_s;
        m_s_tot = (1-eps)*V_MH_tot*rho_s;

        /***************************************************************************************************************************************************************
******************************************************************  SOLID BED TEMPERATURE  *********************************************************************
****************************************************************************************************************************************************************/

         A = ((1-eps)*rho_s*w*cp_g + rho_eff*C_MH)*V_tube;  // Coefficient of Temperature derivative - Considers the effective heat capacity C_MH = eps*rho_g*cp_g + (1-eps)*rho_s*cp_s = 500 J/(kg m3)

         B = ((1-eps) + ((1-eps)*rho_s/rho_g*w+eps)*(1 - a_pg*T_bed) - eps)*V_tube; // Coefficient of pressure derivative

        //  E = ((1-eps)*rho_s * DeltaH/(MW_g/1000))*V_tube;  // Coefficient of gravimetric-density derivative

        // D = (der(m_g_tot_tube)*(h_g_in-h_g) - q_hex); // known term (energy in, and its time variation in the tank, and exchanged heat)
        if Charging == true then
        D = (der(m_g_tot_tube)*(h_g_in - h_g) - q_hex); // known term (energy in, and its time variation in the tank, and exchanged heat)
         E = ((1-eps)*rho_s * DeltaH/(MW_g/1000))*V_tube;  // Coefficient of gravimetric-density derivative

        else
          D = (der(m_g_tot_tube)*(h_g_in - h_g)  - q_hex);    // known term (energy in, and its time variation in the tank, and exchanged heat)
          E = ((1-eps)*rho_s * DeltaH_des/(MW_g/1000))*V_tube;  // Coefficient of gravimetric-density derivative

        end if;

          q_hex = U*A_ext*(T_bed - (Tc_in+Tc_out)/2); // heat exchanged for one tube only
         //   q_hex = U*A_ext*(T_bed - Tc_in); // heat exchanged for one tube only

        //  q_hex = U*A_ext*(DeltaT_lm*Ft); // heat exchanged for one tube only

         U = 1/(((ID/2)/k_eff)*A_in/A_ext + 1/h + R_th_cont +  R_th_tube);

          R_th_tube = OD/2*log(OD/ID)/(k_tube); // [m2.K/W] thermal resistance of the tube wall
        /*
 DeltaT_lm = ((T_bed  - Tc_in) - (T_bed  - Tc_out)) / Modelica.Math.log(DeltaT_lm_DEN);

DeltaT_lm_NUM = ((T_bed - Tc_in)-(T_bed - Tc_out));

 DeltaT_lm_DEN = ((((T_bed - Tc_in)  /(T_bed - Tc_out)))^2)^(1/2);

DeltaT_lm_DEN1 = (T_bed - Tc_in);

DeltaT_lm_DEN2 = (T_bed - Tc_out);
*/
        heatTransferCoefficientPort.h = h;

         Tc_out = q_hex*heatTransferCoefficientPort.Nt/(mDot_c*cp_c) + Tc_in;

         der(T_bed)= (D - B*der(p) - E *der(w))/A; // Energy balance

        /***************************************************************************************************************************************************************
  ***************************************************************  MH properties  ******************************************************************************
  **************************************************************************************************************************************************************/
        rho_eff = (1-eps)*rho_s + eps * rho_g; // effective density

        x_in = F_rp;
        y_out = DATA_interpolation(x=x_in);

        x_k_eff_in = p;
        y_k_eff_out = DATA_interpolation1(x_k_eff=x_k_eff_in);

        if VariableProperties == true then

        C_MH = y_out;
        k_eff = y_k_eff_out;

        else
        k_eff = MH_properties.k_s;
          C_MH = 500;
        end if;

        /*************************************************************************************************************************************************************
***************************************************************  Volume calculations  ***************************************************************************
**************************************************************************************************************************************************************/

        V_MH_tot = heatTransferCoefficientPort.Nt*V_tube;
        V_ratio = V_MH_tot/V_in;
        V_ratio_OD = (pi*OD^2/4*L*heatTransferCoefficientPort.Nt)/V_in;
        /***************************************************************************************************************************************************************
  ***************************************************************  Heat in the tank  ******************************************************************************
  **************************************************************************************************************************************************************/

        q_gen_p_tube = eps*V_tube*der(p); // Heat of compression [1]
        q_gen_DeltaH_tube = -E*der(w); // Heat of reaction
        q_gen_tot_tube = q_gen_p_tube + q_gen_DeltaH_tube;

        q_gen_p_tot = heatTransferCoefficientPort.Nt*eps*V_tube*der(p); // Heat of compression [1]
        q_gen_DeltaH_tot = heatTransferCoefficientPort.Nt*(-E*der(w)); // Heat of reaction
        q_gen_tot = q_gen_p_tot + q_gen_DeltaH_tot;
        /***************************************************************************************************************************************************************
  *******************************************  Thermodynamic states and properties calculations (COOLPROP)  ****************************************************
  **************************************************************************************************************************************************************/

         drhodt = drhodP_T * der(p) + drhodT_P * der(T_bed); // derivative of the hydrogen density

          v_g=1/rho_g;

         hydrogen0=Medium.setState_pT(p=p0, T=Tamb_start);
         hydrogen=Medium.setState_pT(p=p, T=T_bed);
         hydrogen1=Medium.setState_pT(p=p, T=Tamb_start);
         hydrogen2=Medium.setState_pT(p=p+dP, T=T_bed);
         hydrogen3=Medium.setState_pT(p=p, T=T_bed+dT);
         hydrogen4=Medium.setState_pT(p=p_amb, T=Tamb_start);

          h_g_in=Medium.specificEnthalpy(hydrogen0);
         // h_g0=Medium.specificEnthalpy(hydrogen1);

         h_g=Medium.specificEnthalpy(hydrogen) "Dynamic enthalpy";
         cp_g=Medium.specificHeatCapacityCp(hydrogen);
         u_g=Medium.specificInternalEnergy(hydrogen);
          rho_g0=Medium.density(hydrogen0);
            rho_g0_full=Medium.density(hydrogen0);

         rho_g=Medium.density(hydrogen);
         s_g=Medium.specificEntropy(hydrogen);
         cv_g=Medium.specificHeatCapacityCv(hydrogen);
         mu_g=Medium.dynamicViscosity(hydrogen);

         rho_dP =Medium.density(hydrogen2);
         cp_dP= Medium.specificHeatCapacityCp(hydrogen2);

         drhodP_T = (rho_dP - rho_g) / dP;
         drhodT_P = (rho_dT - rho_g) / dT;

         a_pg = rho_g*(1/rho_dT - 1/rho_g) / dT;
         k_tg = - rho_g*(1/rho_dP - 1/rho_g)/dP;

         rho_dT = Medium.density(hydrogen3);
         cv_dT = Medium.specificHeatCapacityCv(hydrogen3);

         dcpdP_T = (cp_dP - cp_g) / dP;
         dcvdT_P = (cv_dT - cv_g) / dT;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************ REFERENCES ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /*
[1] M. Visaria, I. Mudawar, T. Pourpoint, and S. Kumar, “Study of heat transfer and kinetics parameters influencing 
    the design of heat exchangers for hydrogen storage in high-pressure metal hydrides,” Int. J. Heat Mass Transf., vol. 
    53, no. 9–10, pp. 2229–2239, Apr. 2010.
[2] T. L. Pourpoint, V. Velagapudi, I. Mudawar, Y. Zheng, and T. S. Fisher, “Active cooling of a metal hydride 
    system for hydrogen storage,” Int. J. Heat Mass Transf., vol. 53, no. 7–8, pp. 1326–1332, Mar. 2010.
[3]
[4]

*/

         annotation (Diagram(graphics), Icon(graphics={Bitmap(extent={{-64,-50},
                      {66,46}}, fileName=
                      "modelica://HySDeP/Graphics/MH_bed.png")}));
        end MH_bed;
      end Zero_Dimensional_Models;

      annotation (Icon(graphics));
    end Tanks_With_HEX;

    package Ports
      connector p_T

        import SI = Modelica.SIunits;

      SI.Pressure p0;
      SI.Pressure p_amb;
      SI.Temperature T;

        annotation (Icon(graphics={Ellipse(
                extent={{-48,52},{60,-54}},
                lineColor={0,0,0},
                fillColor={0,0,0},
                fillPattern=FillPattern.Solid)}));
      end p_T;

      connector HeatPort

        import SI = Modelica.SIunits;
      type VolumicHeatFlowRate = Real (final quantity="VolumicHeatFlowRate", final unit="W/m3");

      SI.Temperature T;
      flow VolumicHeatFlowRate qDot;

        annotation (Icon(graphics={Ellipse(
                extent={{-54,53},{54,-53}},
                lineColor={0,0,0},
                fillColor={255,0,0},
                fillPattern=FillPattern.Solid)}));

      end HeatPort;

      connector PressurePort

        import SI = Modelica.SIunits;

      SI.Pressure p;

        annotation (Icon(graphics={Ellipse(
                extent={{-48,52},{60,-54}},
                lineColor={0,0,0},
                fillColor={0,127,0},
                fillPattern=FillPattern.Solid)}));
      end PressurePort;

      connector HeatTransferCoefficientPort

        import SI = Modelica.SIunits;

      SI.CoefficientOfHeatTransfer h;
      Integer Nt;
        annotation (Icon(graphics={Ellipse(
                extent={{-48,52},{60,-54}},
                lineColor={0,0,0},
                fillColor={255,255,0},
                fillPattern=FillPattern.Solid), Text(
                extent={{-54,36},{68,-34}},
                lineColor={0,0,0},
                textString="0D")}));
      end HeatTransferCoefficientPort;

      connector HeatTransferCoefficientPort_1D

        import SI = Modelica.SIunits;

      SI.Temperature T;
      flow SI.HeatFlowRate qDot;
      Integer Counter "Counts pieces of walls in the discretized (1D) model";
        annotation (Icon(graphics={Ellipse(
                extent={{-48,52},{60,-54}},
                lineColor={0,0,0},
                fillColor={255,255,0},
                fillPattern=FillPattern.Solid), Text(
                extent={{-38,20},{50,-22}},
                lineColor={0,0,0},
                textString="1D")}));
      end HeatTransferCoefficientPort_1D;

      connector Nt_h

        import SI = Modelica.SIunits;

      Real Nt; // number of tubes
      SI.CoefficientOfHeatTransfer h;

        annotation (Icon(graphics={Ellipse(
                extent={{-48,52},{60,-54}},
                lineColor={0,0,0},
                fillColor={0,0,0},
                fillPattern=FillPattern.Solid), Text(
                extent={{38,-48},{104,-78}},
                lineColor={0,0,0},
                textString="h; Nt")}));
      end Nt_h;

      connector Temperature_port

        import SI = Modelica.SIunits;

      SI.Temperature T;

        annotation (Icon(graphics={Ellipse(
                extent={{-48,52},{60,-54}},
                lineColor={0,0,0},
                fillColor={0,0,0},
                fillPattern=FillPattern.Solid), Text(
                extent={{-24,36},{38,-44}},
                lineColor={0,0,255},
                textString="T")}));
      end Temperature_port;

      connector Nt_Port

        import SI = Modelica.SIunits;

      Integer Nt; // number of tubes

        annotation (Icon(graphics={Ellipse(
                extent={{-48,52},{60,-54}},
                lineColor={0,0,0},
                fillColor={255,170,85},
                fillPattern=FillPattern.Solid), Text(
                extent={{38,-48},{104,-78}},
                lineColor={0,0,0},
                textString="Nt")}));
      end Nt_Port;

      connector Mass_Port

        import SI = Modelica.SIunits;

      SI.Mass m_tot;
      // SI.Mass m_g_tot;
      // SI.Mass m_abs_tot;

        annotation (Icon(graphics={Ellipse(
                extent={{-48,44},{50,-54}},
                lineColor={0,0,0},
                fillColor={170,85,255},
                fillPattern=FillPattern.Solid)}));
      end Mass_Port;
    end Ports;

    package Sources
      model Pressure_ramp_and_Tamb

        import SI = Modelica.SIunits;

      // Call of the Coolprop package to get the H2 properties
       replaceable package Medium =
          ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
          Modelica.Media.Interfaces.PartialMedium;

        /*****************************************************************  PORTS  **********************************************************************/
      Ports.p_T p_T annotation (Placement(transformation(extent={{-8,-72},{12,-52}}),
              iconTransformation(extent={{-10,-66},{10,-46}})));

        Ports.PressurePort p_ramp
          annotation (Placement(transformation(extent={{-98,22},{-78,42}}),
              iconTransformation(extent={{-110,46},{-90,66}})));

        /***************************************************************** VARIABLE DECLARATION **********************************************************************/

        parameter SI.Temperature Tamb = 293.15 "K";
      parameter SI.Pressure  p0 = 101325 "initial pressure in tank [Pa]";
        parameter SI.Pressure height = 3e7 "target pressure in tank"; // same input values of the Visaria and Mudawar paper (Study of ....) for everything: material properties (1.5wt% and thermal and kinetic properties) pressure ramp and alpha(2500 W/m2K))
       parameter SI.Time duration = 60 "Duration of ramp";
      parameter SI.Pressure  p_amb = 101325 "Ambient pressure [Pa]"; // Questa diventerá una inner property insieme a quella del tank (che deve essere uguale!) e verrá definita da qualche altra parte dove questo modello e Tank saranno messi dentro

        parameter SI.Time startTime = 0 "Output = offset for time < startTime";
      SI.Pressure p;

      equation
      p = p0 + (if time < startTime then 0 else if time < startTime + duration then ((time - startTime) * (height - p0)) / duration else height - p0);

      p_ramp.p = p;

       p_T.p0=p0;
       p_T.p_amb=p_amb;
        p_T.T=Tamb;

        annotation (
            Diagram(graphics), Icon(graphics={Bitmap(extent={{-98,-58},{100,72}},
                  fileName=
                    "modelica://HySDeP/Graphics/Hydrogen refueling station - H2 source.jpg")}));
      end Pressure_ramp_and_Tamb;

      model Coolant_Source

        Ports.Temperature_port temperature_port
          annotation (Placement(transformation(extent={{64,-12},{84,8}})));

        import SI = Modelica.SIunits;
      outer parameter SI.Temperature Tc_in = 273.15;

      equation
        temperature_port.T = Tc_in;
        annotation (Icon(graphics={Ellipse(
                extent={{-72,62},{72,-64}},
                lineColor={0,0,255},
                fillColor={170,255,255},
                fillPattern=FillPattern.Solid)}));
      end Coolant_Source;
      annotation (Icon(graphics));
    end Sources;

    package Properties
      package MH
        record BaseProperties

          import SI = Modelica.SIunits;
        type GravimetricDensity = Real(quantity = "Gravimetric Density", unit = "kg/kg");

        /***************** General Properties **************************/
          type MolarEnthalpy = Real(quantity = "MolarEnthalpy", unit = "J/mol");

           parameter GravimetricDensity w_max = 0.015
            "Max gravimetric density for the material [kgH2/kgTOT]";
          parameter MolarEnthalpy DeltaH = -14390
            "Enthalpy of generation [J/mol]";
            parameter MolarEnthalpy DeltaH_des = -24500
            "Enthalpy of generation [J/mol]";                                             // [2]
          parameter SI.MolarEntropy DeltaS = -91.3
            "Entropy of generation [J/(mol*K)]";
            parameter SI.MolarEntropy DeltaS_des = -122
            "Entropy of generation [J/(mol*K)]";
          constant SI.Frequency Ca = 150 "Activation rate [1/s]";
           constant SI.Frequency Ca_des = 300 "Activation rate [1/s]";
          constant MolarEnthalpy Ea = 20700 "Activation Energy [J/molH2]";
            constant MolarEnthalpy Ea_des = 16500 "Activation Energy [J/molH2]";
          parameter SI.Density rho_s = 6200 "Crystalline density [kg/m3]";
          parameter SI.SpecificHeatCapacityAtConstantPressure cp_s0 = 500
            "Specific heat capacity [J/(kg K)]";
          parameter SI.ThermalConductivity k_s = 1
            "Thermal conductivity [W/(m K)]";

          annotation (Icon(graphics={Bitmap(extent={{-82,-54},{80,64}},
                    fileName="modelica://HySDeP/Graphics/MH_properties.jpg")}));

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************ REFERENCES ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /*
[1] M. Visaria, I. Mudawar, T. Pourpoint, and S. Kumar, “Study of heat transfer and kinetics parameters influencing 
    the design of heat exchangers for hydrogen storage in high-pressure metal hydrides,” Int. J. Heat Mass Transf., vol. 
    53, no. 9–10, pp. 2229–2239, Apr. 2010.
    
[2] Milan Visaria and IssamMudawar. Experimental investigation and theoretical modeling of dehydriding process in 
    high-pressure metal hydride hydrogen storage systems. international journal of hydrogen energy, 37(7):5735–5749, 2012.
    */
        end BaseProperties;

        record Ti1_1CrMn

        extends MHSS.Properties.MH.BaseProperties;
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************ REFERENCES ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /*
[1] M. Visaria, I. Mudawar, T. Pourpoint, and S. Kumar, “Study of heat transfer and kinetics parameters influencing 
    the design of heat exchangers for hydrogen storage in high-pressure metal hydrides,” Int. J. Heat Mass Transf., vol. 
    53, no. 9–10, pp. 2229–2239, Apr. 2010.
    */

          annotation (Icon(graphics={Bitmap(extent={{-90,-66},{94,70}},
                    fileName=
                      "modelica://HySDeP/Graphics/Ti11CrMn_properties.jpg")}));
        end Ti1_1CrMn;

        record Ti1_1CrMn_with_Enhanced_Conductivity

        extends MHSS.Properties.MH.BaseProperties(
            w_max=0.015 "Max gravimetric density for the material [kgH2/kgTOT]",
            DeltaH=-14390 "Enthalpy of generation [J/mol]",
            DeltaS=-91.3 "Entropy of generation [J/(mol*K)]",
            Ca=150 "Activation rate [1/s]",
            Ea=20700 "Activation Energy [J/molH2]",
            rho_s=6200 "Crystalline density [kg/m3]",
            cp_s0=500 "Specific heat capacity [J/(kg K)]",
            k_s=10 "Thermal conductivity [W/(m K)]",
            DeltaH_des=-24500 "Enthalpy of generation [J/mol]",
            DeltaS_des=-112 "Entropy of generation [J/(mol*K)]",
            Ca_des=300 "Activation rate [1/s]",
            Ea_des=16500 "Activation Energy [J/molH2]");

          annotation (Icon(graphics={Bitmap(extent={{-78,-60},{80,70}},
                    fileName=
                      "modelica://HySDeP/Graphics/Ti11CrMn_enhanced_conductivity_properties.jpg")}));
        end Ti1_1CrMn_with_Enhanced_Conductivity;

        record LaNi5

        extends MHSS.Properties.MH.BaseProperties(
            w_max=0.015 "Max gravimetric density for the material [kgH2/kgTOT]",
            DeltaH=-14390 "Enthalpy of generation [J/mol]",
            DeltaS=-91.3 "Entropy of generation [J/(mol*K)]",
            Ca=150 "Activation rate [1/s]",
            Ea=20700 "Activation Energy [J/molH2]",
            rho_s=6200 "Crystalline density [kg/m3]",
            cp_s0=500 "Specific heat capacity [J/(kg K)]",
            k_s=1 "Thermal conductivity [W/(m K)]",
            DeltaH_des=-24500 "Enthalpy of generation [J/mol]",
            DeltaS_des=-112 "Entropy of generation [J/(mol*K)]",
            Ca_des=300 "Activation rate [1/s]",
            Ea_des=16500 "Activation Energy [J/molH2]");

          annotation (Icon(graphics={Bitmap(extent={{-64,-64},{72,66}},
                    fileName="modelica://HySDeP/Graphics/Lani5_properties.jpg")}));
        end LaNi5;
      end MH;

      package Coolant
        partial model Coolant
          // Partial model including cooling fluid  properties, to be extended by filling mmodels and by HEAT TRANSFER COEFFICIENT calculating model
          // in another partial model:
          import SI = Modelica.SIunits;
          /* *************************************************************************************************
****************A TABLE FOR PROPERTIES AT DIFFERENT TEMPERATURES (INITIAL CONDITIONS) **************
*******************************************MUST BE DONE********************************************* */
          // Here referring to Dexcool (all data from "SYSTEM SIMULATION MODEL FOR HIGH PRESSURE METAL HYDRIDE HYDROGEN STORAGE SYSTEMS" paper)
          // ***********************  T=273 K  *********************
          parameter SI.Temperature T_c = 273 "K";
          // two temperatures: T=273 and T= 338
          parameter SI.Density rho_c = 1060 "kg/m3";
          parameter SI.SpecificHeatCapacity c_c = 3460 "J/(kg K)";
          parameter SI.ThermalConductivity k_c = 0.415 "W/(m K)";
          parameter SI.KinematicViscosity mu_c = 6.8e-6 "m2/s";
          parameter SI.VolumeFlowRate Q_c = 40 "L/min/tube";
          /*
// ***********************  T=338 K  *********************

parameter SI.Temperature T_c=338 "K";   // two temperatures: T=273 and T= 338

parameter SI.Density rho_c=1037 "kg/m3";
parameter SI.SpecificHeatCapacity c_c= 3730   "J/(kg K)";
parameter SI.ThermalConductivity k_c=0.459 "W/(m K)";
parameter SI.KinematicViscosity m_c = 1.1+e-06 "m2/s";
parameter SI.VolumeFlowRate m_c=40 "L/min/tube";

*/
        end Coolant;

        record CoolantBaseProperties

          import SI = Modelica.SIunits;

        /***************** General Properties of DexCool at 273.15 K [1] **************************/

        parameter SI.Density rho_c = 1068 "kg/m3";
        parameter SI.SpecificHeatCapacity cp_c = 3460 "J/(kg K)";
        parameter SI.ThermalConductivity k_c = 0.415 "W/(m K)";
        parameter SI.DynamicViscosity mu_c = 7.26e-3 "Pa s";

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************ REFERENCES ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /*
[1] M. Visaria, I. Mudawar, and T. Pourpoint, “Enhanced heat exchanger design for hydrogen 
storage using high-pressure metal hydride: Part 1. Design methodology and computational 
results,” Int. J. Heat Mass Transf., vol. 54, no. 1–3, pp. 413–423, Jan. 2011.
    */
          annotation (Icon(graphics={Bitmap(extent={{-82,-66},{84,76}},
                    fileName=
                      "modelica://HySDeP/Graphics/CoolantBaseProperties.png")}));
        end CoolantBaseProperties;

        record Dex_Cool

        extends CoolantBaseProperties;
          import SI = Modelica.SIunits;

        /***************** General Properties of DexCool at 273.15 K [1] **************************/
        /*
parameter SI.Density rho_c = 1068 "kg/m3";
parameter SI.SpecificHeatCapacity cp_c = 3460 "J/(kg K)";
parameter SI.ThermalConductivity k_c = 0.415 "W/(m K)";
parameter SI.DynamicViscosity mu_c = 7.26e-3 "Pa s";
*/

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************ REFERENCES ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /*
[1] M. Visaria, I. Mudawar, and T. Pourpoint, “Enhanced heat exchanger design for hydrogen 
storage using high-pressure metal hydride: Part 1. Design methodology and computational 
results,” Int. J. Heat Mass Transf., vol. 54, no. 1–3, pp. 413–423, Jan. 2011.
    */

          annotation (Icon(graphics={Bitmap(extent={{-90,-76},{88,86}},
                    fileName="modelica://HySDeP/Graphics/Dex-Cool.png")}));
        end Dex_Cool;

        record Kerosene_example_Book_by_Serth

          import SI = Modelica.SIunits;
        extends CoolantBaseProperties(
        rho_c = 784.905 "kg/m3",
        cp_c = 2470.2 "J/(kg K)",
        k_c = 0.137 "W/(m K)",
        mu_c = 0.00040098 "Pa s");
        /***************** General Properties of DexCool at 273.15 K [1] **************************/

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************ REFERENCES ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /*
[1] M. Visaria, I. Mudawar, and T. Pourpoint, “Enhanced heat exchanger design for hydrogen 
storage using high-pressure metal hydride: Part 1. Design methodology and computational 
results,” Int. J. Heat Mass Transf., vol. 54, no. 1–3, pp. 413–423, Jan. 2011.
    */

        end Kerosene_example_Book_by_Serth;
      end Coolant;

      package Tube
        record BaseProperties

          import SI = Modelica.SIunits;

        /***************** General Properties **************************/

          parameter SI.ThermalConductivity k_tube = 187
            "Thermal conductivity [W/(m K)]";                                             // [1]
          parameter SI.Density rho_tube = 2700 "Thermal conductivity [W/(m K)]";          // [1]
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************ REFERENCES ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /*
[1] http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA6061t6
    */
        end BaseProperties;

        record Aluminum6061

        extends MHSS.Properties.Tube.BaseProperties(k_tube=187
              "Thermal conductivity [W/(m K)]", rho_tube = 2700
              "Density [kg/m3]");
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************ REFERENCES ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /*
[1] http://www.engineeringtoolbox.com/heat-exchanger-material-thermal-conductivities-d_1488.html
    
    ISO/TR 15916:2004 Basic considerations for the safety of hydrogen systems /Technical Report: (most of) Aluminum alloys; Coper alloys and austenitic steels do not present embrittlement (or this phenomenon is very modest)
    
    */

          annotation (Icon(graphics={Bitmap(extent={{-76,-70},{90,74}},
                    fileName=
                      "modelica://HySDeP/Graphics/Table Aluminum 6061.jpg")}));
        end Aluminum6061;

        record AISI_304_stainless_steel

        extends MHSS.Properties.Tube.BaseProperties(k_tube=16
              "Thermal conductivity [W/(m K)]", rho_tube = 8000
              "Density [kg/m3]");
                                                           // [1]

          annotation (Icon(graphics={Bitmap(extent={{-68,-66},{78,78}},
                    fileName="modelica://HySDeP/Graphics/Table AISI 304.jpg")}));
        end AISI_304_stainless_steel;

        record Aluminum

        extends MHSS.Properties.Tube.BaseProperties(k_tube=202
              "Thermal conductivity [W/(m K)]", rho_tube = 2700
              "Density [kg/m3]");
                                                           // [1]
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************ REFERENCES ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /*
[1] http://www.engineeringtoolbox.com/heat-exchanger-material-thermal-conductivities-d_1488.html
    
    ISO/TR 15916:2004 Basic considerations for the safety of hydrogen systems /Technical Report: (most of) Aluminum alloys; Coper alloys and austenitic steels do not present embrittlement (or this phenomenon is very modest)
    
    */

          annotation (Icon(graphics={Bitmap(extent={{-76,-70},{90,74}},
                    fileName="modelica://HySDeP/Graphics/Table Aluminum.jpg")}));
        end Aluminum;

        record Copper

        extends MHSS.Properties.Tube.BaseProperties(k_tube=386
              "Thermal conductivity [W/(m K)]", rho_tube = 8940
              "Density [kg/m3]");

          annotation (Icon(graphics={Bitmap(extent={{-70,-76},{88,74}},
                    fileName="modelica://HySDeP/Graphics/Table Copper.jpg")}));
        end Copper;

      end Tube;
    end Properties;

    package Heat_Transfer

      package Heat_Exchangers
        model Tube_In_Tube

          import SI = Modelica.SIunits;
          import Modelica.Constants.pi;

         outer parameter Boolean Calculated_Volume
            "if true then the volume is calculated from L_input and d_input";

        outer parameter SI.Length delta;
        outer parameter Real a_r "Aspect ratio; 0<a<1 (ID/OD)"; // outer
        outer parameter SI.MassFlowRate mDot_c;  // mass flow rate in the HEX (in the annulus)
        SI.MassFlowRate mDot_c_tube;  // mass flow rate per tube (in the annulus)

        final parameter SI.Length ID = 2*delta "inner tube's inner diameter";
        outer parameter Real PR;
        final parameter SI.Length OD_tube_in_tube = ID/a_r;  // outer tube's inner diameter
        final parameter SI.Length Dh = OD_tube_in_tube-ID;  // hydraulic diameter
        final parameter SI.Length CTP=0.93;  // tube count constant for one tube pass [1]
        Real CL;
        outer parameter Modelica.SIunits.Conversions.NonSIunits.Angle_deg theta_tp
            "Layout angle (30; 45; 60; 90)";
        outer parameter SI.Length d_tank
            "tank diameter - used if Calculated_Volume=true";
        outer parameter SI.Volume V_input
            "Inner volume of cylindrical tank - used if Calculated_Volume=false";

        SI.Length Ds = (if Calculated_Volume == false then (4*V_input/(Modelica.Constants.pi*Ltube))^(1/2) else d_tank)
            "tank diameter";
        outer parameter SI.Length Ltube = 1 "Tube active length"; // outer = Ltube
        final parameter SI.Length ODo = OD_tube_in_tube/0.8; // Outer tube's outer diameter (for Nt_max_tube_in_tube calculation) => considers tube's thickness
        final parameter Real fg = 1.615*(1 + 0.14*a_r^(-0.5));

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  PROPERTIES  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        outer parameter SI.Density rho_c;
        outer parameter SI.SpecificHeatCapacity cp_c;
        outer parameter SI.ThermalConductivity k_c;
        outer parameter SI.DynamicViscosity mu_c;
        SI.KinematicViscosity nu_c;

        Integer Nt_max_tube_in_tube;

        SI.Area A_cross_tube_in_tube; // Cross flow area
        Real Re_tube_in_tube;
        Real Pr;
        SI.Velocity v_c_tube_in_tube;

        // VARIABLES FOR Nu CALCULATION IN LAMINAR REGIME
        Real Nu1;  // Nu for thermally developed flow
        Real Nu11;  // Nu for thermally developing flow
        Real Nu111;

        // VARIABLES FOR Nu CALCULATION IN TRANSITION REGIME
        Real Nu11_Lam_Tran; // Nu for thermally developing flow
        Real Nu111_Lam_Tran;
        Real Nu_Lam_Tran; // Mean Nu for hydrodinamic developing and thermally developing flow

        Real k_Tur_Tran;  // factor for turb. regime
        Real Re1_Tur_Tran;
        Real f_Tur_Tran; // friction factor for turbulent regime
        Real Nu_Tur_Tran;
        Real g;

        //  VARIABLES FOR Nu CALCULATION IN TURBULENT REGIME
        Real k;  // factor for turbulent regime
        Real Re_1;
        Real f;  // friction factor (turbulent regime)
        Real f_corr;     // Annular Tube Nu correction Factor

        Real Nu;
        SI.CoefficientOfHeatTransfer h_tube_in_tube;

        equation
        A_cross_tube_in_tube = (pi/4*OD_tube_in_tube^2) - (pi/4*ID^2); // [m2] Cross flow area
        v_c_tube_in_tube = mDot_c_tube/(A_cross_tube_in_tube*rho_c); // coolant velocity
        nu_c = mu_c/rho_c;
        Pr =rho_c*nu_c*cp_c/k_c;     // Prandtl number
        Re_tube_in_tube = v_c_tube_in_tube*Dh/nu_c;           // Reynolds number

        if (theta_tp == 30 or theta_tp ==60) then

         CL = 0.87;
        else  // it refers to the case of theta_tp == 45deg or 90deg
         CL = 1;
        end if;

        mDot_c_tube = mDot_c/Nt_max_tube_in_tube; // mass flow rate per tube (in the annulus)
        Nt_max_tube_in_tube = integer(floor(0.785 * (CTP/CL)* Ds^2/(PR^2*ODo^2) + 0.5)); // Eq. 8.9, Ref. [1]

        // Laminar
        Nu1 = 3.66 + 1.2*a_r^(-0.8); // Nu for thermally developed flow
        Nu11 =  fg*(Re_tube_in_tube*Pr*Dh/Ltube)^(1/3); // Nu for thermally developing flow
        Nu111 = (2/(1 + 22*Pr))^(1/6)*(Re_tube_in_tube*Pr*Dh/Ltube)^(1/2);

        // Transition
        Nu11_Lam_Tran =  fg*(2300*Pr*Dh/Ltube)^(1/3); // Nu for thermally developing flow
        Nu111_Lam_Tran = (2/(1 + 22*Pr))^(1/6)*(2300*Pr*Dh/Ltube)^(1/2);
        Nu_Lam_Tran = (Nu1^3 + Nu11_Lam_Tran^3 + Nu111_Lam_Tran^3)^(1/3); // Mean Nu for hydrodinamic developing and thermally developing flow
        k_Tur_Tran = 1.07 + 900/10e3 - 0.63/(1 + 10*Pr); // factor for turb. regime
        Re1_Tur_Tran = 10e3* ((1+a_r^2)*log(a_r) + (1 - a_r^2))/((1-a_r)^2*log(a_r));
        f_Tur_Tran = (1.8* log10(Re1_Tur_Tran) - 1.5)^(-2); // friction factor for turbulent regime
        Nu_Tur_Tran = f_corr*((f_Tur_Tran/8)*10e3*Pr)/(k_Tur_Tran + 12.7*(f_Tur_Tran/8)^(0.5)
                                                                                            *(Pr^(2/3) - 1))*(1+(Dh/Ltube)^(2/3));
        g = (Re_tube_in_tube - 2300)/(1e4 - 2300);

        // Turbulent
        k = 1.07 + 900/Re_tube_in_tube - 0.63/(1 + 10*Pr); // factor for turbulent regime
        Re_1 = Re_tube_in_tube* ((1+a_r^2)*log(a_r) + (1 - a_r^2))/((1-a_r)^2*log(a_r));
        f = (1.8 * log10(Re_1) - 1.5)^(-2); // friction factor (turbulent regime)
        f_corr = 0.75*a_r^(-0.17);     // Annular Tube Nu correction Factor

        if Re_tube_in_tube<2300 then
        Nu = (Nu1^3 + Nu11^3 + Nu111^3)^(1/3); // Mean Nu for hydrodinamic developing and thermally developing flow
        else if (Re_tube_in_tube>=2300 and Re_tube_in_tube<1e4) then
        Nu = (1-g)*Nu_Lam_Tran + g*Nu_Tur_Tran;

        else
        Nu = f_corr*((f/8)*Re_tube_in_tube*Pr)/(k + 12.7*(f/8)^(0.5)
                                                                   *(Pr^(2/3) - 1))*(1+(Dh/Ltube)^(2/3)); // Modified Gnielinski, Ref. [5]
        end if;
          end if;

        h_tube_in_tube = k_c*Nu/Dh; // heat convection coefficient

          annotation (Icon(graphics={Bitmap(extent={{-62,-56},{70,64}},
                    fileName=
                      "modelica://HySDeP/Graphics/tube-in-tube-heat-exchanger.jpg")}));
        end Tube_In_Tube;

        model Shell_and_tube

          import SI = Modelica.SIunits;
          import Modelica.Constants.pi;

        type ShellSideMassVelocity = Real(quantity= "ShellSideMassVelocity", unit="kg/(s.m2)");

        /******************************  GRAPHICS  ************************************ */

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  BOOLEAN PARAMETERS  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        parameter Boolean NumberOfTubes_Given = false;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  HEAT TRANSFER PARAMETERS  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        parameter SI.MassFlowRate mDot_c = 0.3 "Coolant mass flow rate";   // OUTER
        parameter Real Ft = 0.875 "temperature factor (1 = countercurrent)"; // minimum value is 0.75 outer

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  GEOMETRY PARAMETERS  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        // PARAMETERS
        parameter SI.Length ID = 2*delta "inner tube diameter = 2*delta";
        parameter SI.Length Ltube = 1 "Tube length";    // active tube length THIS MUST BE EQUAL TO THE TANK LENGTH // OUTER
        parameter Real PR = 1.25 "Tube pitch to tube diameter ratio: Pt/OD";  // Tube pitch to tube outer diameter ratio (Pt/OD) http://www.engineeringpage.com/technology/thermal/pitch.html
        parameter SI.Length Ds = 0.4 "Shell Diameter";    // THIS MUST BE EQUAL TO THE TANK INNER DIAMETER // OUTER
        parameter Modelica.SIunits.Conversions.NonSIunits.Angle_deg theta_tp = 30
            "Layout angle (30,60,45,90deg)";                                                                        // in deg
        parameter Integer Nt_input = 100 "Used if NumberOfTubes_Given = true";
        parameter Real Nss = 1
            "Number of Sealing strips (pairs) in one baffle spacing";
        parameter Real B0=0.5
            "Initial value for the baffle space 0<B<1 - percentage (>0.2 (>5cm in general))";                       // Baffle space as percentage of the shell diameter "Ds"  // OUTER
        parameter Real Bc = 25e-2 "Baffle cut (20-49%) with 20-25% optimum";
        parameter Integer Nb = 4
            "Number of baffles - input if NumberOfTubes_Given = true";
        parameter Real Lo_star = 1
            "1=<Lo_star<2 --> Ratio between the inlet baffle spacing and the central baffle spacing"; // Ratio between the outlet baffle spacing and the central baffle spacing
        parameter Real Li_star = 1
            "1=<Li_star<2 --> Ratio between the inlet baffle spacing and the central baffle spacing"; // Ratio between the inlet baffle spacing and the central baffle spacing
        parameter Integer N_shell = 1; // Number of shells
        //outer parameter SI.Temperature T_in = 273.15 "Inlet coolant temperature";

        // FINAL PARAMETERS
        final parameter SI.Length OD = ID/0.8 "inner tube diameter";
        final parameter SI.Length delta = 1e-2 "inner tube diameter";
        final parameter SI.Length Pt = PR*OD; // tube pitch
        final parameter SI.Length Cl_shell=0.025;   /* [m] - Distance between shell and tube bundle. From: "Mechanical Constraints of thermal design of Shell and Tube Exchangers" Available at: https://www.academia.edu */
        final parameter SI.Length CTP=0.93;  // tube count constant for one tube pass [1]
        final parameter SI.Length Cl_tube = Pt - OD; // Distance between tubes
        final parameter SI.Length Ratio = Cl_tube./Pt;
        final parameter SI.Length B_space0= B0*Ds; // [m] - Baffle space in absolute value to be used initially. the actual one is recalculate = B_space
        final parameter Real Nb_real= Ltube/B_space0 - 1;  // number of baffles (real number = decimal)
        final parameter SI.Length Lpl=0;  //  =0 for single tube pass. Expresses the effect of the tube lane partition bypass width;
        final parameter SI.Length Lbb = (12*1e-3 + 0.005*Ds); // [m] Bundle to shell clearance
        final parameter SI.Length D_otl = Ds-Lbb; // Tube bank OTL diameter
        final parameter SI.Length D_ctl = D_otl - OD; // Bundle diameter
        final parameter SI.Length Ltb = 0.4e-3; // Diametral clearance between tube outside diameter and baffle hole
        final parameter SI.Length Lsb = (3.1*1e-3 + 0.004*Ds); // [7]

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  Constants and Constant exponents  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        constant Real a1_1 = 1.4;
        constant Real a1_2 = 1.36;
        constant Real a1_3 = 0.593;
        constant Real a1_4 = 0.321;
        constant Real a2_1 = -0.667;
        constant Real a2_2 = -0.657;
        constant Real a2_3 = -0.477;
        constant Real a2_4 = -0.388;
        constant Real g_c = 1; // [6]

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  VARIABLES **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        // Geometry
        SI.Length Dn;  // Nozzles diameter from [6] Ch. 5 (See table 5.3)
        Real x;
        Real r_s;
        Real r_lm;
        Real f_s; // Friction factor [1]
        SI.Length Lts; // Tube sheet thickness
        SI.Length Lto; // Shell length
        SI.Length B_space; // Actual baffle space (recalculated)
        SI.Length Ltp_eff;
        SI.Area A_shell_BD;
        SI.Angle theta_ds; // in rad to be used in the formulae
        SI.Angle theta_ctl; // in rad to be used in the formulae
        Real Fw;  // Fraction of tubes in the baffle window
        Real Fc;  // Fraction of tubes in pure cross flow
        SI.Area Swg; // Gross window flow area
        SI.Area Swt; // Segmental baffle window area
        SI.Area Sw; // Net crossflow area throughout one baffle window
        SI.Length Dw;  // Equivalent hydraulic diameter of a segmental baffle window (used only for DeltaP calculations)
        SI.Length Lpp; // Effective tube row distance in the flow direction (function of theta_tp)
        Integer Nb_integer;  // number of baffles (rounded number = integer)
        Real Re_s; // Reynolds number (shell side)
        Real CL;   /* tube layout constant (value valid for inclination angle by 30deg and 60deg) =1 for 45 or 90deg => more tube density for 30 and 60deg. See Ref. [1] */
        Real B;     // Baffle space as percentage of the shell diameter "Ds"   (B=B_space/Ds)
        Real a;
        Real a1;
        Real a2;
        Real a3;
        Real a4;
        Real b1;
        Real b2;
        Real b3;
        Real b4;
        Real b;
        SI.Velocity vf_Kern; // for water and similar fluids it should be in the range: 0.6-1.5 m/s "Heat Exchanger Design Book"
        SI.Length De; // equivalent diameter
        SI.CoefficientOfHeatTransfer h_BD;
        SI.CoefficientOfHeatTransfer h_kern;
        SI.CoefficientOfHeatTransfer h_id;
        SI.Area A_shell_kern;
        Integer Nt_max; // maximum number of tubes that fit in the shell
        Integer Nt_w; // Number of tubes in window area
        Real Nt_cc; // Number of tube rows in crossflow (i.e. between tbaffle tips)
        Real Nt_cw; // Effective number of tube rows crossed in the baffle window
        Integer Nc; //  total number of tube rows crossed in the entire heat exchanger
        SI.Area Ahex;
        SI.Area Sb; // Bypass area between the shell and the tube bundle within one baffle
        SI.Area Stb; // Tube to baffle hole leakagearea for one baffle (to calculate baffle leakage effect parameters J1 and R1)
        SI.Area Ssb; // Shell-to-baffle leakage area for one baffle (to calculate baffle leakage effect parameters J1 and R1)
        Real Fsbp; // Ratio of the bypass area (Sb) to the overall crossflow area (Sm) - used for J1 and R1
        ShellSideMassVelocity Gs_kern; // Shell side mass velocity
        ShellSideMassVelocity Gs_BD; // Shell side mass velocity
        ShellSideMassVelocity Gw; // Window mass velocity
        ShellSideMassVelocity Gn; // Nozzles mass velocity
        Real Nu_kern;
        Real r_ss; // =Nss/Nt_cc
        Real C_bh;
        Real C_bp;
        Real n_s; // exponent in "Js" (0.6 ia a value valid for turbulent flow)
        Real n_b; // exponent in "Js" (0.6 ia a value valid for turbulent flow)

        // Properties
        Real Re_n; // Reynolds number in the nozzles
        Real Pr; // Prandtl number

        // Correction factors (J -> heat transfer; R -> Pressure losses)
        Real Ji; // Colburn factor
        Real J1; // Correction factor for baffle leakage effects, including both shell-to-baffle and tube-to-baffle leakage.
        Real Jc; // Correction factor for segmental baffle window (considers baffle cut and baffle spacing)
        Real Jr; // Only applicable for Laminar flow (Re_s<=100 and fullt effective for Re_s<20) with the maximum limit Jr=0.4
        Real Jb; // Correction factors for bundle bypass effects for heat transfer
        Real Js; // Heat transfer correction factor for unequal baffle spacing at inlet and/or outlet
        Real J_tot; // Cumulative correction factor (Good design: J_tot>=0.6)

        Real Rb; // Correction factors for bundle bypass effects for pressure drop
        Real R1;
        Real Rs; // Pressure drop correction factor for unequal baffle spacing at inlet and/or outlet

        // Pressure Losses
        SI.Pressure DeltaP_id; // Ideal shell-side tube-bank pressure drop
        SI.Pressure DeltaP_tot_no_nozzles; // Total pressure dorp on the shell side neglecting the nozzles
        SI.Pressure DeltaP_tot_with_nozzles; // Total pressure dorp on the shell side including the nozzles
        SI.Pressure DeltaP_c; // Pressure drop in all central baffle spaces
        SI.Pressure DeltaP_w; // Pressure drop in one window
        SI.Pressure DeltaP_e; // Pressure drop in the entrance and exit baffle spaces
        SI.Pressure DeltaP_w_tot; // Pressure drop in all windowst baffle spaces
        SI.Pressure DeltaP_n; // Pressure drop in the shell's nozzles

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  PROPERTIES  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        inner replaceable parameter
            MHSS.Properties.Coolant.CoolantBaseProperties Coolant
            constrainedby MHSS.Properties.Coolant.CoolantBaseProperties                                         annotation (choicesAllMatching=true);

        final parameter SI.Density rho_c = Coolant.rho_c "kg/m3";
        final parameter SI.SpecificHeatCapacity cp_c = Coolant.cp_c "J/(kg K)";
        final parameter SI.ThermalConductivity k_c = Coolant.k_c "W/(m K)";
        final parameter SI.DynamicViscosity mu_c = Coolant.mu_c "Pa s";

        constant SI.Density rho_w = 1000;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  EQUATIONS  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        equation
        if NumberOfTubes_Given == false then
          Nt_max = integer(floor(0.785 * (CTP/CL)* Ds^2/(PR^2*OD^2) + 0.5)); // Eq. 8.9, Ref. [1]
          Nb_integer = if (Nb_real > 0) then integer(floor(Nb_real + 0.5)) else integer(ceil(Nb_real - 0.5));
        else
          Nt_max = Nt_input;
            Nb_integer = Nb;

        end if;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  KERN METHOD  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        A_shell_kern = (B0*Ds^2*Cl_tube)/Pt; /* [m2] Cross flow area shell-side; 
  modified equation (original Eq 8.15 pp. 309 [1]: 
  Ds*Cl_tube*Distance_Baffle/Pt, here rearranged with this substitution: 
  Distance_Baffle = B*Ds where: 20%< B <100% ). Note that Cl_tube/Pt =Ratio
  that is: Cl_tube/Pt = (Pt-OD)/Pt = 1-1/PR and thus is a constant.
  At the end for each simulation where Ds and B are constants, A_shell_kern is a
  constant as well. */

        vf_Kern = mDot_c/(rho_c*A_shell_kern);
        Ahex = Modelica.Constants.pi*Ltube*OD*Nt_max;
        Pr = (cp_c*mu_c)/k_c;
        Gs_kern = mDot_c/A_shell_kern;
        De = 4*(Pt^2 - Modelica.Constants.pi*OD^2/4)/(Modelica.Constants.pi*OD); // shell hydraulic (or equivalent) diameter - square layout
        Nu_kern = 0.36*(De*Gs_kern/mu_c)^0.55 * (cp_c*mu_c/k_c)^(1/3);
        h_kern = k_c/De*Nu_kern;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  BELL-DELAWARE METHOD **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        /**************************************************  GEOMETRY CALCULATIONS  ********************************************************************************************************************************************************************************************************************************************/

        A_shell_BD = B_space*(Lbb + D_ctl/Ltp_eff*(Pt-OD));
        Lto = Ltube + 2*Lts;
        B_space = Ltube/(Nb_integer + 1);  // exact baffle spacing
        B = B_space/Ds;
        h_id = Ji*cp_c*Gs_BD*(1/Pr)^(2/3); //  the parameter that considers the different dynamic viscosity at the wall and in the bulk fluid according to the temperature is not considered here yet (compute it when you have the temperature)
        Re_s = OD*mDot_c/(mu_c*A_shell_BD); // used for the Bell-Delaware equation, Ref. [1]
        Gs_BD = mDot_c/A_shell_BD;

        theta_ds = 2*(acos(1-2*Bc)); // Here "B" is given already in decimal form, tehrefore do not write B/100 as in the original formula of the book
        theta_ctl = 2*(acos(Ds/D_ctl*(1-2*Bc)));

        Fw = theta_ctl/(2*pi) - sin(theta_ctl)/(2*pi);
        Fc = 1-2*Fw;

        Swg = pi/4*Ds^2*(theta_ds/(2*pi) - sin(theta_ds)/(2*pi));
        Swt = Nt_max*Fw*pi/4*OD^2;
        Sw = Swg-Swt;
        Sb = B_space*(Ds - D_otl + Lpl);
        Ssb = pi/2*Ds*Lsb*((2*pi-theta_ds)/(2*pi));
        Stb = pi/4*((OD + Ltb)^2 - OD^2)*Nt_max*(1-Fw); // [7, 8]

        Nt_w = integer(floor(Nt_max*Fw + 0.5));
        // Nt_cc = integer(floor(Ds/Lpp*(1-2*Bc) + 0.5));
        Nt_cc = Ds/Lpp*(1-2*Bc);
        // Nt_cw = integer(floor(0.8/Lpp*(Ds*Bc - (Ds - D_ctl)/2) + 0.5));
        Nt_cw = 0.8/Lpp*(Ds*Bc - (Ds - D_ctl)/2);
        Nc = integer(floor((Nt_cc+Nt_cw)*(Nb_integer+1) + 0.5));

        Fsbp= Sb/A_shell_BD;
        Dw = 4*Sw/(pi*OD*Nt_w + pi*Ds*theta_ds/(2*pi));

        if theta_tp == 30 then
        Lpp = 0.866*Pt;
        else if  (theta_tp == 60) then
         Lpp = 0.866*Pt; // TO BE CHECKED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             else if (theta_tp == 90) then
           Lpp = Pt;
               else
        Lpp = 0.707*Pt; //  value for theta_tp = 45deg
        end if;
        end if;
        end if;

        if Ds*0.1>25e-3 then
        Lts = 0.1*Ds;
        else
          Lts = 25.4e-3;
        end if;

        if (theta_tp == 30 or theta_tp ==90) then
        Ltp_eff = Pt;
        else if  (theta_tp == 60) then
         Ltp_eff = 0.866*Pt; // TO BE CHECKED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          else
        Ltp_eff = 0.707*Pt; //  value for theta_tp = 45deg
        end if;
        end if;

        x = (-0.15*(1+r_s) + 0.8);

        r_s = Ssb/(Ssb + Stb);

        r_lm = (Ssb + Stb)/A_shell_BD;

        r_ss = Nss/Nt_cc;

        f_s = b1*(1.33/(Pt/OD))^b*(Re_s)^(b2); // [1]

        Gn = mDot_c/(pi*Dn^2/4);

        Re_n = 4*mDot_c/(pi*Dn*mu_c);

        Gw = mDot_c/(A_shell_BD*Sw)^(1/2);

        /**************************************************  COEFFICIENTS AND EXPONENTS  ********************************************************************************************************************************************************************************************************************************************/

        if (theta_tp == 30 or theta_tp ==60) then
          a3 = 1.450;
          a4 = 0.519;
          b3 = 7.00;
          b4 = 0.5;
         CL = 0.87;
        if Re_s<10 then
        a1 = a1_1;
        a2 = a2_1;
        b1= 48;
        b2 = -1;
        else if Re_s >= 10 and Re_s < 100 then
        a1= a1_2;
        a2= a2_2;
        b1= 45.1;
        b2 = -0.973;

        else if Re_s>=100 and Re_s<1000 then
        a1 = a1_3;
        a2 = a2_3;
        b1= 4.570;
        b2 = -0.476;

             else if Re_s>=1000 and Re_s<10000 then
        a1 = a1_4;
        a2 = a2_4;
        b1= 0.486;
        b2 = -0.152;

             else
        a1 = a1_4;
        a2 = a2_4;
        b1= 0.372;
        b2 = -0.123;

        end if;
        end if;
        end if;
        end if;
        else if theta_tp == 45 then
          a3 = 1.930;
          a4 = 0.5;
          b3 = 6.59;
          b4 = 0.520;
        CL = 1;

            if Re_s<10 then
        a1 = 1.550;
        a2 = -0.667;
        b1= 32;
        b2 = -1;

        else if Re_s >= 10 and Re_s < 100 then
        a1= 0.498;
        a2= -0.656;
        b1= 26.2;
        b2 = -0.913;

        else if Re_s>=100 and Re_s<1000 then
        a1 = 0.730;
        a2 = -0.5;
        b1= 3.5;
        b2 = -0.476;

             else  if Re_s>=1000 and Re_s<10000 then
        a1 = 0.37;
        a2 = -0.396;
        b1= 0.333;
        b2 = -0.136;

             else
            a1 = 0.37;
        a2 = -0.396;
        b1= 0.303;
        b2 = -0.126;

        end if;
        end if;
        end if;
        end if;

        else    // so this is for the case theta_tp = 90deg
          a3 = 1.187;
          a4 = 0.370;
          b3 = 6.30;
          b4 = 0.378;
        CL = 1;

          if Re_s<10 then
        a1 = 0.970;
        a2 = -0.667;
        b1= 35;
        b2 = -1;

        else if Re_s >= 10 and Re_s < 100 then
        a1= 0.90;
        a2= -0.631;
        b1= 32.1;
        b2 = -0.963;

        else if Re_s>=100 and Re_s<1000 then
        a1 = 0.408;
        a2 = -0.46;
        b1= 6.09;
        b2 = -0.602;
             else if Re_s>=1000 and Re_s<10000 then
        a1 = 0.107;
        a2 = -0.266;
        b1= 0.0815;
        b2 = +0.022;

             else
        a1 = 0.370;
        a2 = -0.395;
        b1= 0.391;
        b2 = -0.148;

             end if;
        end if;
          end if;
          end if;
        end if;
        end if;

        a = a3/(1 + 0.14*(Re_s)^(a4));
        b = b3/(1 + 0.14*(Re_s)^b4);

        if r_ss< 0.5 then
        Rb = exp(-C_bp*Fsbp*(1-(2*r_ss)^(1/3)));
        else
        Rb = min(exp(-C_bp*Fsbp*(1-(2*r_ss)^(1/3))), 1);
        end if;

        if Re_s >= 100 then
        C_bh = 1.35;
        C_bp = 3.7;
        else
        C_bh = 1.25;
        C_bp = 4.5;
        end if;

        if Ds>=4*0.0254 and Ds<=10*0.0254 then
          Dn = 2*0.0254;
        elseif Ds>=12*0.0254 and Ds<=17.5*0.0254 then
          Dn = 3*0.0254;
          elseif Ds>=19.25*0.0254 and Ds<=21.25*0.0254 then
          Dn = 4*0.0254;
          elseif Ds>=23*0.0254 and Ds<=29*0.0254 then
          Dn = 6*0.0254;

          elseif Ds>=31*0.0254 and Ds<=37*0.0254 then
          Dn = 8*0.0254;

          elseif Ds>=39*0.0254 and Ds<=42*0.0254 then
          Dn = 10*0.0254;
        else
          Dn = 12*0.0254;

        end if;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  Heat Transfer and Pressure-drop correction factors  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        Ji = a1*(1.33/PR)^a*Re_s^(a2); // Colburn coefficient
        Jc = 0.55 + 0.72*Fc;
        J1 = 0.44*(1-r_s) + (1-0.44*(1-r_s))*Modelica.Constants.e^(-2.2*r_lm);
        Jb = min(exp(-C_bh*Fsbp*(1 - (2*r_ss)^(1/3))), 1);
        J_tot = Jc*J1*Jb*Js*Jr;

        if Re_s >= 100 then
        Jr = 1;
          n_s = 0.6;
          n_b = 0.2;
        Js = ((Nb_integer - 1) + Li_star^(1-n_s) + Lo_star^(1-n_s))/((Nb_integer - 1) + (Li_star) + (Lo_star)); // [6] by Serth
        DeltaP_w = (2+0.6*Nt_cw)*Gw^2/(2*g_c*rho_c);

        elseif (20<Re_s and Re_s<100) then
          Jr = max( 1.51/(Nc^(0.18)) + ((20-Re_s)/80)*(1.51/(Nc^(0.18)) - 1), 0.4); // maximum limit for Jr
           n_s = 1/3;
          n_b = 1;
        Js = (((Nb_integer - 1) + Li_star^(1-n_s) + Lo_star^(1-n_s))/((Nb_integer - 1) + (Li_star) + (Lo_star)) + 1) /2; // Foor laminar flow Js is about halfwas between 1 and the value calculated for turbulent conditions
        DeltaP_w = 26*Gw*mu_c/(g_c*rho_c)*(Nt_cw/(Pt-OD) + B_space/(Dw^2)) + 2*Gw^2/(g_c*rho_c);
        else
            Jr = max( 1.51/(Nc^(0.18)), 0.4); // maximum limit for Jr
             n_s = 1/3;
          n_b = 1;
        Js = (((Nb_integer - 1) + Li_star^(1-n_s) + Lo_star^(1-n_s))/((Nb_integer - 1) + (Li_star) + (Lo_star)) + 1) /2; // Foor laminar flow Js is about halfwas between 1 and the value calculated for turbulent conditions
        DeltaP_w = 26*Gw*mu_c/(g_c*rho_c)*(Nt_cw/(Pt-OD) + B_space/(Dw^2)) + 2*Gw^2/(g_c*rho_c);
        end if;

        R1 = Modelica.Constants.e^(-1.33*(1+r_s)*r_lm^(x));
        Rs =  ((1/Li_star)^(2-n_b) + (1/Lo_star)^(2-n_b)); // [7, 8]

        /***************************************************************  Heat Transfer coefficient (Bell_Delaware)  *******************************************************************************************************************************************************************************************************************************/
        h_BD = h_id*Jc*J1*Jb*Js*Jr;

         /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /***************************************************************  PRESSURE DROPS  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        DeltaP_id = (2*f_s*Nt_cc*(Gs_BD)^(2))/(g_c*rho_c); // Ideal tube bank pressure drop
        DeltaP_c = (Nb_integer -1)*(DeltaP_id*Rb*R1); // Combined pressure drop of all the interior crossflow sections
        DeltaP_e = 2*DeltaP_id*(1 + Nt_cw/Nt_cc)*Rb*Rs; // Combined pressure drop for the entrance and exit sections
        // DeltaP_w = pressure loss in windows depends on Re_s, thus is inserted in an "if-statement"
        DeltaP_w_tot = DeltaP_w*Nb_integer*R1; // Pressure drop in all windows
        DeltaP_tot_no_nozzles = DeltaP_e+ DeltaP_w_tot + DeltaP_c; // Total Pressure drop (neglecting pressure drop in the nozzles)

          if Re_n >= 100 then
            DeltaP_n = 2e-13*N_shell*Gn^2/(rho_c/rho_w);

        else
            DeltaP_n = 4e-13*N_shell*Gn^2/(rho_c/rho_w);
        end if;

        DeltaP_tot_with_nozzles = DeltaP_tot_no_nozzles + DeltaP_n;
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************ REFERENCES ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /*
[1] S. Kakac, H. Liu. "Heat Exchangers Selection, rating and Thermal 
    Design". 2nd Edition, CRC Press. ISBN: 0-8493-0902-6.
  
[2] Serth R.W. Process Heat Transfer, Principles and Applications. Ed. 2007 ISBN: 978-0-12-373588-1

[3] Heat Exchanger Design Book

[4] Wolverine Tube Heat Transfer Data Book
    
*/

          annotation (Icon(graphics={Bitmap(extent={{-90,-56},{94,78}},
                    fileName=
                      "modelica://HySDeP/Graphics/shell_and_tube_heat_exchanger_flow.png")}));
        end Shell_and_tube;

        record Heat_Exchanger_Configuration

          import SI = Modelica.SIunits;

        /***************** Record of the parameter "HEX" used to select the heat exchanger configuration **************************/

        parameter Integer HEX = 1;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************ REFERENCES ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /*
[1] M. Visaria, I. Mudawar, and T. Pourpoint, “Enhanced heat exchanger design for hydrogen 
storage using high-pressure metal hydride: Part 1. Design methodology and computational 
results,” Int. J. Heat Mass Transf., vol. 54, no. 1–3, pp. 413–423, Jan. 2011.
    */

          annotation (Icon(graphics={Bitmap(extent={{-82,-86},{74,90}},
                    fileName=
                      "modelica://HySDeP/Graphics/Record_Heat_Exchangers_LIST.png")}));
        end Heat_Exchanger_Configuration;

        record SHELL_AND_TUBE

          import SI = Modelica.SIunits;

        /***************** Record of the parameter "HEX" used to select the heat exchanger configuration **************************/

        extends Heat_Exchanger_Configuration(
        HEX = 1);
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************ REFERENCES ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /*
[1] M. Visaria, I. Mudawar, and T. Pourpoint, “Enhanced heat exchanger design for hydrogen 
storage using high-pressure metal hydride: Part 1. Design methodology and computational 
results,” Int. J. Heat Mass Transf., vol. 54, no. 1–3, pp. 413–423, Jan. 2011.
    */

        end SHELL_AND_TUBE;

        record TUBE_IN_TUBE

          import SI = Modelica.SIunits;

        /***************** Record of the parameter "HEX" used to select the heat exchanger configuration **************************/
        extends Heat_Exchanger_Configuration(
        HEX = 2);

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************ REFERENCES ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /*
[1] M. Visaria, I. Mudawar, and T. Pourpoint, “Enhanced heat exchanger design for hydrogen 
storage using high-pressure metal hydride: Part 1. Design methodology and computational 
results,” Int. J. Heat Mass Transf., vol. 54, no. 1–3, pp. 413–423, Jan. 2011.
    */

        end TUBE_IN_TUBE;

        model Shell_and_tube_With_Port

          import SI = Modelica.SIunits;
          import Modelica.Constants.pi;

        type ShellSideMassVelocity = Real(quantity= "ShellSideMassVelocity", unit="kg/(s.m2)");

        /******************************  GRAPHICS  ************************************ */

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  BOOLEAN PARAMETERS  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        parameter Boolean NumberOfTubes_Given = false;
        outer parameter Boolean Bell_Delaware_Method = true;
        outer parameter Boolean Calculated_Volume
            "if true then the volume is calculated from L_input and d_input";
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  HEAT TRANSFER PARAMETERS  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        parameter SI.MassFlowRate mDot_c = 30 "Coolant mass flow rate";   // OUTER
        parameter Real Ft = 0.875 "temperature factor (1 = countercurrent)"; // minimum value is 0.75

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  GEOMETRY PARAMETERS  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        // PARAMETERS (Inner in "Tank System")

        outer parameter SI.Length d_tank
            "tank diameter - used if Calculated_Volume=true";
        outer parameter SI.Volume V_input
            "Inner volume of cylindrical tank - used if Calculated_Volume=false";

        SI.Length Ds = (if Calculated_Volume == false then (4*V_input/(Modelica.Constants.pi*Ltube))^(1/2) else d_tank)
            "tank diameter";

        outer parameter SI.Length ID "inner tube diameter = 2*delta";
        outer parameter SI.Length Ltube = 1 "Tube length";
        outer parameter Real PR = 1.25
            "Tube pitch to tube diameter ratio: Pt/OD";                             // Tube pitch to tube outer diameter ratio (Pt/OD) http://www.engineeringpage.com/technology/thermal/pitch.html
        // outer parameter SI.Length Ds = 0.4 "Shell Diameter";
        outer parameter Modelica.SIunits.Conversions.NonSIunits.Angle_deg theta_tp = 30
            "Layout angle (30,60,45,90deg)";                                                                        // in deg
        outer parameter Integer Nt_input = 100
            "Used if NumberOfTubes_Given = true";
        outer parameter Real Nss = 1
            "Number of Sealing strips (pairs) in one baffle spacing";
        outer parameter Real B0=0.5
            "Initial value for the baffle space 0<B<1 - percentage (>0.2 (>5cm in general))";                       // Baffle space as percentage of the shell diameter "Ds"  // OUTER
        outer parameter Real Bc = 25e-2
            "Baffle cut (20-49%) with 20-25% optimum";
        outer parameter Integer Nb = 4
            "Number of baffles - input if NumberOfTubes_Given = true";
        outer parameter Real Lo_star = 1
            "1=<Lo_star<2 --> Ratio between the inlet baffle spacing and the central baffle spacing"; // Ratio between the outlet baffle spacing and the central baffle spacing
        outer parameter Real Li_star = 1
            "1=<Li_star<2 --> Ratio between the inlet baffle spacing and the central baffle spacing"; // Ratio between the inlet baffle spacing and the central baffle spacing
        outer parameter Integer N_shell = 1; // Number of shells

        // FINAL PARAMETERS
        final parameter SI.Length OD = ID/0.8 "inner tube diameter";
        final parameter SI.Length Pt = PR*OD; // tube pitch
        final parameter SI.Length Cl_shell=0.025;   /* [m] - Distance between shell and tube bundle. From: "Mechanical Constraints of thermal design of Shell and Tube Exchangers" Available at: https://www.academia.edu */
        final parameter SI.Length CTP=0.93;  // tube count constant for one tube pass [1]
        final parameter SI.Length Cl_tube = Pt - OD; // Distance between tubes
        final parameter SI.Length Ratio = Cl_tube./Pt;
        final parameter SI.Length Lpl=0;  //  =0 for single tube pass. Expresses the effect of the tube lane partition bypass width;
        final parameter SI.Length Ltb = 0.4e-3; // Diametral clearance between tube outside diameter and baffle hole

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  Constants and Constant exponents  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        constant Real a1_1 = 1.4;
        constant Real a1_2 = 1.36;
        constant Real a1_3 = 0.593;
        constant Real a1_4 = 0.321;
        constant Real a2_1 = -0.667;
        constant Real a2_2 = -0.657;
        constant Real a2_3 = -0.477;
        constant Real a2_4 = -0.388;
        constant Real g_c = 1; // [6]

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  VARIABLES **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        // Geometry
        SI.Length Dn;  // Nozzles diameter from [6] Ch. 5 (See table 5.3)
        Real x;
        Real r_s;
        Real r_lm;
        Real f_s; // Friction factor [1]
        SI.Length Lts; // Tube sheet thickness
        SI.Length Lto; // Shell length
        SI.Length B_space; // Actual baffle space (recalculated)
        SI.Length Ltp_eff;
        SI.Area A_shell_BD;
        SI.Angle theta_ds; // in rad to be used in the formulae
        SI.Angle theta_ctl; // in rad to be used in the formulae
        Real Fw;  // Fraction of tubes in the baffle window
        Real Fc;  // Fraction of tubes in pure cross flow
        SI.Area Swg; // Gross window flow area
        SI.Area Swt; // Segmental baffle window area
        SI.Area Sw; // Net crossflow area throughout one baffle window
        SI.Length Dw;  // Equivalent hydraulic diameter of a segmental baffle window (used only for DeltaP calculations)
        SI.Length Lpp; // Effective tube row distance in the flow direction (function of theta_tp)
        Integer Nb_integer;  // number of baffles (rounded number = integer)
        Real Re_s; // Reynolds number (shell side)
        Real CL;   /* tube layout constant (value valid for inclination angle by 30deg and 60deg) =1 for 45 or 90deg => more tube density for 30 and 60deg. See Ref. [1] */
        Real B;     // Baffle space as percentage of the shell diameter "Ds"   (B=B_space/Ds)
        Real a;
        Real a1;
        Real a2;
        Real a3;
        Real a4;
        Real b1;
        Real b2;
        Real b3;
        Real b4;
        Real b;
        SI.Velocity vf_Kern; // for water and similar fluids it should be in the range: 0.6-1.5 m/s "Heat Exchanger Design Book"
        SI.Length De; // equivalent diameter
        SI.CoefficientOfHeatTransfer h_BD;
        SI.CoefficientOfHeatTransfer h_kern;
        SI.CoefficientOfHeatTransfer h_id;
        SI.CoefficientOfHeatTransfer h = (if Bell_Delaware_Method == false then h_kern else h_BD);

        SI.Area A_shell_kern;
        Integer Nt_max; // maximum number of tubes that fit in the shell
        Integer Nt_w; // Number of tubes in window area
        Real Nt_cc; // Number of tube rows in crossflow (i.e. between tbaffle tips)
        Real Nt_cw; // Effective number of tube rows crossed in the baffle window
        Integer Nc; //  total number of tube rows crossed in the entire heat exchanger
        SI.Area Ahex;
        SI.Area Sb; // Bypass area between the shell and the tube bundle within one baffle
        SI.Area Stb; // Tube to baffle hole leakagearea for one baffle (to calculate baffle leakage effect parameters J1 and R1)
        SI.Area Ssb; // Shell-to-baffle leakage area for one baffle (to calculate baffle leakage effect parameters J1 and R1)
        Real Fsbp; // Ratio of the bypass area (Sb) to the overall crossflow area (Sm) - used for J1 and R1
        ShellSideMassVelocity Gs_kern; // Shell side mass velocity
        ShellSideMassVelocity Gs_BD; // Shell side mass velocity
        ShellSideMassVelocity Gw; // Window mass velocity
        ShellSideMassVelocity Gn; // Nozzles mass velocity
        Real Nu_kern;
        Real r_ss; // =Nss/Nt_cc
        Real C_bh;
        Real C_bp;
        Real n_s; // exponent in "Js" (0.6 ia a value valid for turbulent flow)
        Real n_b; // exponent in "Js" (0.6 ia a value valid for turbulent flow)
        SI.Length B_space0; // [m] - Baffle space in absolute value to be used initially. the actual one is recalculate = B_space
        SI.Length Lbb;  // [m] Bundle to shell clearance
        SI.Length D_otl; // Tube bank OTL diameter
        SI.Length Lsb; // [7]
        SI.Length D_ctl; // Bundle diameter
        Real Nb_real;   // number of baffles (real number = decimal)

        // Properties
        Real Re_n; // Reynolds number in the nozzles
        Real Pr; // Prandtl number

        // Correction factors (J -> heat transfer; R -> Pressure losses)
        Real Ji; // Colburn factor
        Real J1; // Correction factor for baffle leakage effects, including both shell-to-baffle and tube-to-baffle leakage.
        Real Jc; // Correction factor for segmental baffle window (considers baffle cut and baffle spacing)
        Real Jr; // Only applicable for Laminar flow (Re_s<=100 and fullt effective for Re_s<20) with the maximum limit Jr=0.4
        Real Jb; // Correction factors for bundle bypass effects for heat transfer
        Real Js; // Heat transfer correction factor for unequal baffle spacing at inlet and/or outlet
        Real J_tot; // Cumulative correction factor (Good design: J_tot>=0.6)

        Real Rb; // Correction factors for bundle bypass effects for pressure drop
        Real R1;
        Real Rs; // Pressure drop correction factor for unequal baffle spacing at inlet and/or outlet

        // Pressure Losses
        SI.Pressure DeltaP_id; // Ideal shell-side tube-bank pressure drop
        SI.Pressure DeltaP_tot_no_nozzles; // Total pressure dorp on the shell side neglecting the nozzles
        SI.Pressure DeltaP_tot_with_nozzles; // Total pressure dorp on the shell side including the nozzles
        SI.Pressure DeltaP_c; // Pressure drop in all central baffle spaces
        SI.Pressure DeltaP_w; // Pressure drop in one window
        SI.Pressure DeltaP_e; // Pressure drop in the entrance and exit baffle spaces
        SI.Pressure DeltaP_w_tot; // Pressure drop in all windowst baffle spaces
        SI.Pressure DeltaP_n; // Pressure drop in the shell's nozzles

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  PROPERTIES  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        outer parameter SI.Density rho_c;
        outer parameter SI.SpecificHeatCapacity cp_c;
        outer parameter SI.ThermalConductivity k_c;
        outer parameter SI.DynamicViscosity mu_c;

        constant SI.Density rho_w = 1000;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  HEAT EXCHANGER TYPE  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        // Metal hydride properties stored in records [1]:
        inner replaceable parameter
            MHSS.Heat_Transfer.Heat_Exchangers.Heat_Exchanger_Configuration HEX
            constrainedby
            MHSS.Heat_Transfer.Heat_Exchangers.Heat_Exchanger_Configuration                                     annotation (choicesAllMatching=true);

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  EQUATIONS  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

          Ports.HeatTransferCoefficientPort heatTransferCoefficientPort annotation (Placement(transformation(extent={{-6,52},{14,72}}), iconTransformation(extent={{-14,44},{14,72}})));

        equation
        if NumberOfTubes_Given == false then
          Nt_max = integer(floor(0.785 * (CTP/CL)* Ds^2/(PR^2*OD^2) + 0.5)); // Eq. 8.9, Ref. [1]
          Nb_integer = if (Nb_real > 0) then integer(floor(Nb_real + 0.5)) else integer(ceil(Nb_real - 0.5));
        else
          Nt_max = Nt_input;
            Nb_integer = Nb;

        end if;

         heatTransferCoefficientPort.h = h;
          heatTransferCoefficientPort.Nt = Nt_max;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  KERN METHOD  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        A_shell_kern = (B0*Ds^2*Cl_tube)/Pt; /* [m2] Cross flow area shell-side; 
  modified equation (original Eq 8.15 pp. 309 [1]: 
  Ds*Cl_tube*Distance_Baffle/Pt, here rearranged with this substitution: 
  Distance_Baffle = B*Ds where: 20%< B <100% ). Note that Cl_tube/Pt =Ratio
  that is: Cl_tube/Pt = (Pt-OD)/Pt = 1-1/PR and thus is a constant.
  At the end for each simulation where Ds and B are constants, A_shell_kern is a
  constant as well. */

        vf_Kern = mDot_c/(rho_c*A_shell_kern);
        Ahex = Modelica.Constants.pi*Ltube*OD*Nt_max;
        Pr = (cp_c*mu_c)/k_c;
        Gs_kern = mDot_c/A_shell_kern;
        De = 4*(Pt^2 - Modelica.Constants.pi*OD^2/4)/(Modelica.Constants.pi*OD); // shell hydraulic (or equivalent) diameter - square layout
        Nu_kern = 0.36*(De*Gs_kern/mu_c)^0.55 * (cp_c*mu_c/k_c)^(1/3);
        h_kern = k_c/De*Nu_kern;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  BELL-DELAWARE METHOD **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        /**************************************************  GEOMETRY CALCULATIONS  ********************************************************************************************************************************************************************************************************************************************/

        B_space0= B0*Ds;
        Lbb = (12*1e-3 + 0.005*Ds); // [m] Bundle to shell clearance
        D_otl = Ds-Lbb; // Tube bank OTL diameter
        Lsb = (3.1*1e-3 + 0.004*Ds); // [7]
        Nb_real= Ltube/B_space0 - 1;  // number of baffles (real number = decimal)
        D_ctl = D_otl - OD; // Bundle diameter

        A_shell_BD = B_space*(Lbb + D_ctl/Ltp_eff*(Pt-OD));
        Lto = Ltube + 2*Lts;
        B_space = Ltube/(Nb_integer + 1);  // exact baffle spacing
        B = B_space/Ds;
        h_id = Ji*cp_c*Gs_BD*(1/Pr)^(2/3); //  the parameter that considers the different dynamic viscosity at the wall and in the bulk fluid according to the temperature is not considered here yet (compute it when you have the temperature)
        Re_s = OD*mDot_c/(mu_c*A_shell_BD); // used for the Bell-Delaware equation, Ref. [1]
        Gs_BD = mDot_c/A_shell_BD;

        theta_ds = 2*(acos(1-2*Bc)); // Here "B" is given already in decimal form, tehrefore do not write B/100 as in the original formula of the book
        theta_ctl = 2*(acos(Ds/D_ctl*(1-2*Bc)));

        Fw = theta_ctl/(2*pi) - sin(theta_ctl)/(2*pi);
        Fc = 1-2*Fw;

        Swg = pi/4*Ds^2*(theta_ds/(2*pi) - sin(theta_ds)/(2*pi));
        Swt = Nt_max*Fw*pi/4*OD^2;
        Sw = Swg-Swt;
        Sb = B_space*(Ds - D_otl + Lpl);
        Ssb = pi/2*Ds*Lsb*((2*pi-theta_ds)/(2*pi));
        Stb = pi/4*((OD + Ltb)^2 - OD^2)*Nt_max*(1-Fw); // [7, 8]

        Nt_w = integer(floor(Nt_max*Fw + 0.5));
        // Nt_cc = integer(floor(Ds/Lpp*(1-2*Bc) + 0.5));
        Nt_cc = Ds/Lpp*(1-2*Bc);
        // Nt_cw = integer(floor(0.8/Lpp*(Ds*Bc - (Ds - D_ctl)/2) + 0.5));
        Nt_cw = 0.8/Lpp*(Ds*Bc - (Ds - D_ctl)/2);
        Nc = integer(floor((Nt_cc+Nt_cw)*(Nb_integer+1) + 0.5));

        Fsbp= Sb/A_shell_BD;
        Dw = 4*Sw/(pi*OD*Nt_w + pi*Ds*theta_ds/(2*pi));

        if theta_tp == 30 then
        Lpp = 0.866*Pt;
        else if  (theta_tp == 60) then
         Lpp = 0.866*Pt; // TO BE CHECKED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             else if (theta_tp == 90) then
           Lpp = Pt;
               else
        Lpp = 0.707*Pt; //  value for theta_tp = 45deg
        end if;
        end if;
        end if;

        if Ds*0.1>25e-3 then
        Lts = 0.1*Ds;
        else
          Lts = 25.4e-3;
        end if;

        if (theta_tp == 30 or theta_tp ==90) then
        Ltp_eff = Pt;
        else if  (theta_tp == 60) then
         Ltp_eff = 0.866*Pt;
          else
        Ltp_eff = 0.707*Pt; //  value for theta_tp = 45deg
        end if;
        end if;

        x = (-0.15*(1+r_s) + 0.8);

        r_s = Ssb/(Ssb + Stb);

        r_lm = (Ssb + Stb)/A_shell_BD;

        r_ss = Nss/Nt_cc;

        f_s = b1*(1.33/(Pt/OD))^b*(Re_s)^(b2); // [1]

        Gn = mDot_c/(pi*Dn^2/4);

        Re_n = 4*mDot_c/(pi*Dn*mu_c);

        Gw = mDot_c/(A_shell_BD*Sw)^(1/2);

        /**************************************************  COEFFICIENTS AND EXPONENTS  ********************************************************************************************************************************************************************************************************************************************/

        if (theta_tp == 30 or theta_tp ==60) then
          a3 = 1.450;
          a4 = 0.519;
          b3 = 7.00;
          b4 = 0.5;
         CL = 0.87;
        if Re_s<10 then
        a1 = a1_1;
        a2 = a2_1;
        b1= 48;
        b2 = -1;
        else if Re_s >= 10 and Re_s < 100 then
        a1= a1_2;
        a2= a2_2;
        b1= 45.1;
        b2 = -0.973;

        else if Re_s>=100 and Re_s<1000 then
        a1 = a1_3;
        a2 = a2_3;
        b1= 4.570;
        b2 = -0.476;

             else if Re_s>=1000 and Re_s<10000 then
        a1 = a1_4;
        a2 = a2_4;
        b1= 0.486;
        b2 = -0.152;

             else
        a1 = a1_4;
        a2 = a2_4;
        b1= 0.372;
        b2 = -0.123;

        end if;
        end if;
        end if;
        end if;
        else if theta_tp == 45 then
          a3 = 1.930;
          a4 = 0.5;
          b3 = 6.59;
          b4 = 0.520;
        CL = 1;

            if Re_s<10 then
        a1 = 1.550;
        a2 = -0.667;
        b1= 32;
        b2 = -1;

        else if Re_s >= 10 and Re_s < 100 then
        a1= 0.498;
        a2= -0.656;
        b1= 26.2;
        b2 = -0.913;

        else if Re_s>=100 and Re_s<1000 then
        a1 = 0.730;
        a2 = -0.5;
        b1= 3.5;
        b2 = -0.476;

             else  if Re_s>=1000 and Re_s<10000 then
        a1 = 0.37;
        a2 = -0.396;
        b1= 0.333;
        b2 = -0.136;

             else
            a1 = 0.37;
        a2 = -0.396;
        b1= 0.303;
        b2 = -0.126;

        end if;
        end if;
        end if;
        end if;

        else    // so this is for the case theta_tp = 90deg
          a3 = 1.187;
          a4 = 0.370;
          b3 = 6.30;
          b4 = 0.378;
        CL = 1;

          if Re_s<10 then
        a1 = 0.970;
        a2 = -0.667;
        b1= 35;
        b2 = -1;

        else if Re_s >= 10 and Re_s < 100 then
        a1= 0.90;
        a2= -0.631;
        b1= 32.1;
        b2 = -0.963;

        else if Re_s>=100 and Re_s<1000 then
        a1 = 0.408;
        a2 = -0.46;
        b1= 6.09;
        b2 = -0.602;
             else if Re_s>=1000 and Re_s<10000 then
        a1 = 0.107;
        a2 = -0.266;
        b1= 0.0815;
        b2 = +0.022;

             else
        a1 = 0.370;
        a2 = -0.395;
        b1= 0.391;
        b2 = -0.148;

             end if;
        end if;
          end if;
          end if;
        end if;
        end if;

        a = a3/(1 + 0.14*(Re_s)^(a4));
        b = b3/(1 + 0.14*(Re_s)^b4);

        if r_ss< 0.5 then
        Rb = exp(-C_bp*Fsbp*(1-(2*r_ss)^(1/3)));
        else
        Rb = min(exp(-C_bp*Fsbp*(1-(2*r_ss)^(1/3))), 1);
        end if;

        if Re_s >= 100 then
        C_bh = 1.35;
        C_bp = 3.7;
        else
        C_bh = 1.25;
        C_bp = 4.5;
        end if;

        if Ds>=4*0.0254 and Ds<=10*0.0254 then
          Dn = 2*0.0254;
        elseif Ds>=12*0.0254 and Ds<=17.5*0.0254 then
          Dn = 3*0.0254;
          elseif Ds>=19.25*0.0254 and Ds<=21.25*0.0254 then
          Dn = 4*0.0254;
          elseif Ds>=23*0.0254 and Ds<=29*0.0254 then
          Dn = 6*0.0254;

          elseif Ds>=31*0.0254 and Ds<=37*0.0254 then
          Dn = 8*0.0254;

          elseif Ds>=39*0.0254 and Ds<=42*0.0254 then
          Dn = 10*0.0254;
        else
          Dn = 12*0.0254;

        end if;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  Heat Transfer and Pressure-drop correction factors  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        Ji = a1*(1.33/PR)^a*Re_s^(a2); // Colburn coefficient
        Jc = 0.55 + 0.72*Fc;
        J1 = 0.44*(1-r_s) + (1-0.44*(1-r_s))*Modelica.Constants.e^(-2.2*r_lm);
        Jb = min(exp(-C_bh*Fsbp*(1 - (2*r_ss)^(1/3))), 1);
        J_tot = Jc*J1*Jb*Js*Jr;

        if Re_s >= 100 then
        Jr = 1;
          n_s = 0.6;
          n_b = 0.2;
        Js = ((Nb_integer - 1) + Li_star^(1-n_s) + Lo_star^(1-n_s))/((Nb_integer - 1) + (Li_star) + (Lo_star)); // [6] by Serth
        DeltaP_w = (2+0.6*Nt_cw)*Gw^2/(2*g_c*rho_c);

        elseif (20<Re_s and Re_s<100) then
          Jr = max( 1.51/(Nc^(0.18)) + ((20-Re_s)/80)*(1.51/(Nc^(0.18)) - 1), 0.4); // maximum limit for Jr
           n_s = 1/3;
          n_b = 1;
        Js = (((Nb_integer - 1) + Li_star^(1-n_s) + Lo_star^(1-n_s))/((Nb_integer - 1) + (Li_star) + (Lo_star)) + 1) /2; // Foor laminar flow Js is about halfwas between 1 and the value calculated for turbulent conditions
        DeltaP_w = 26*Gw*mu_c/(g_c*rho_c)*(Nt_cw/(Pt-OD) + B_space/(Dw^2)) + 2*Gw^2/(g_c*rho_c);
        else
            Jr = max( 1.51/(Nc^(0.18)), 0.4); // maximum limit for Jr
             n_s = 1/3;
          n_b = 1;
        Js = (((Nb_integer - 1) + Li_star^(1-n_s) + Lo_star^(1-n_s))/((Nb_integer - 1) + (Li_star) + (Lo_star)) + 1) /2; // Foor laminar flow Js is about halfwas between 1 and the value calculated for turbulent conditions
        DeltaP_w = 26*Gw*mu_c/(g_c*rho_c)*(Nt_cw/(Pt-OD) + B_space/(Dw^2)) + 2*Gw^2/(g_c*rho_c);
        end if;

        R1 = Modelica.Constants.e^(-1.33*(1+r_s)*r_lm^(x));
        Rs =  ((1/Li_star)^(2-n_b) + (1/Lo_star)^(2-n_b)); // [7, 8]

        /***************************************************************  Heat Transfer coefficient (Bell_Delaware)  *******************************************************************************************************************************************************************************************************************************/
        h_BD = h_id*Jc*J1*Jb*Js*Jr;

         /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /***************************************************************  PRESSURE DROPS  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        DeltaP_id = (2*f_s*Nt_cc*(Gs_BD)^(2))/(g_c*rho_c); // Ideal tube bank pressure drop
        DeltaP_c = (Nb_integer -1)*(DeltaP_id*Rb*R1); // Combined pressure drop of all the interior crossflow sections
        DeltaP_e = 2*DeltaP_id*(1 + Nt_cw/Nt_cc)*Rb*Rs; // Combined pressure drop for the entrance and exit sections
        // DeltaP_w = pressure loss in windows depends on Re_s, thus is inserted in an "if-statement"
        DeltaP_w_tot = DeltaP_w*Nb_integer*R1; // Pressure drop in all windows
        DeltaP_tot_no_nozzles = DeltaP_e+ DeltaP_w_tot + DeltaP_c; // Total Pressure drop (neglecting pressure drop in the nozzles)

          if Re_n >= 100 then
            DeltaP_n = 2e-13*N_shell*Gn^2/(rho_c/rho_w);

        else
            DeltaP_n = 4e-13*N_shell*Gn^2/(rho_c/rho_w);
        end if;

        DeltaP_tot_with_nozzles = DeltaP_tot_no_nozzles + DeltaP_n;
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************ REFERENCES ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /*
[1] S. Kakac, H. Liu. "Heat Exchangers Selection, rating and Thermal 
    Design". 2nd Edition, CRC Press. ISBN: 0-8493-0902-6.
  
[2] Serth R.W. Process Heat Transfer, Principles and Applications. Ed. 2007 ISBN: 978-0-12-373588-1

[3] Heat Exchanger Design Book

[4] Wolverine Tube Heat Transfer Data Book
    
*/

          annotation (Icon(graphics={Bitmap(extent={{-90,-56},{94,78}},
                    fileName=
                      "modelica://HySDeP/Graphics/shell_and_tube_heat_exchanger_flow.png")}));
        end Shell_and_tube_With_Port;

        model Tube_In_Tube_to_validate_with_VDIAtalas

          import SI = Modelica.SIunits;
          import Modelica.Constants.pi;

         outer parameter Boolean Calculated_Volume
            "if true then the volume is calculated from L_input and d_input";

        /** qui devo metter i boolean per il volume, calculated or not ecc ***/
        parameter SI.Length delta = 0.5e-2;
        parameter Real a_r = 0.5 "Aspect ratio; 0<a<1 (ID/OD)"; // outer
        parameter SI.MassFlowRate mDot_c_tot = 30; // outer
        SI.MassFlowRate mDot_c;  // mass flow rate per tube (in the annulus)

        final parameter SI.Length ID = 2*delta "inner tube's inner diameter";
        outer parameter Real PR = 1.25;
        final parameter SI.Length OD = ID/a_r;  // outer tube's inner diameter
        final parameter SI.Length Dh = OD-ID;  // hydraulic diameter
        final parameter SI.Length CTP=0.93;  // tube count constant for one tube pass [1]
        Real CL;
        outer parameter Modelica.SIunits.Conversions.NonSIunits.Angle_deg theta_tp = 30
            "Layout angle (30; 45; 60; 90)";
        outer parameter SI.Length d_tank
            "tank diameter - used if Calculated_Volume=true";
        outer parameter SI.Volume V_input
            "Inner volume of cylindrical tank - used if Calculated_Volume=false";

        SI.Length Ds = (if Calculated_Volume == false then (4*V_input/(Modelica.Constants.pi*Ltube))^(1/2) else d_tank)
            "tank diameter";
        outer parameter SI.Length Ltube = 1 "Tube active length"; // outer = Ltube
        final parameter SI.Length ODo = OD/0.8; // Outer tube's outer diameter (for Nt_max calculation) => considers tube's thickness
        final parameter Real fg = 1.615*(1 + 0.14*a_r^(-0.5));

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  PROPERTIES  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        outer parameter SI.Density rho_c;
        outer parameter SI.SpecificHeatCapacity cp_c;
        outer parameter SI.ThermalConductivity k_c;
        outer parameter SI.DynamicViscosity mu_c;
        SI.KinematicViscosity nu_c;

        Integer Nt_max;

        // Integer Nt ; // number of tubes in integer (i.e. floor...) // ha senso che questo sia l'Nt max calcolabile
        // clearances
        // layout angles etc...

        SI.Area A_cross; // Cross flow area
        Real Re;
        Real Pr;
        SI.Velocity v_c;

        // VARIABLES FOR Nu CALCULATION IN LAMINAR REGIME
        Real Nu1;  // Nu for thermally developed flow
        Real Nu11;  // Nu for thermally developing flow
        Real Nu111;

        // VARIABLES FOR Nu CALCULATION IN TRANSITION REGIME
        Real Nu11_Lam_Tran; // Nu for thermally developing flow
        Real Nu111_Lam_Tran;
        Real Nu_Lam_Tran; // Mean Nu for hydrodinamic developing and thermally developing flow

        Real k_Tur_Tran;  // factor for turb. regime
        Real Re1_Tur_Tran;
        Real f_Tur_Tran; // friction factor for turbulent regime
        Real Nu_Tur_Tran;
        Real g;

        //  VARIABLES FOR Nu CALCULATION IN TURBULENT REGIME
        Real k;  // factor for turbulent regime
        Real Re_1;
        Real f;  // friction factor (turbulent regime)
        Real f_corr;     // Annular Tube Nu correction Factor

        Real Nu;
        SI.CoefficientOfHeatTransfer h;

        equation
        A_cross = (pi/4*OD^2) - (pi/4*ID^2); // [m2] Cross flow area
        v_c = mDot_c/(A_cross*rho_c); // coolant velocity
        nu_c = mu_c/rho_c;
        Pr =rho_c*nu_c*cp_c/k_c;     // Prandtl number
        Re = v_c*Dh/nu_c;           // Reynolds number

        if (theta_tp == 30 or theta_tp ==60) then

         CL = 0.87;
        else  // it refers to the case of theta_tp == 45deg or 90deg
         CL = 1;
        end if;

        mDot_c = mDot_c_tot/Nt_max; // mass flow rate per tube (in the annulus)
        Nt_max = integer(floor(0.785 * (CTP/CL)* Ds^2/(PR^2*ODo^2) + 0.5)); // Eq. 8.9, Ref. [1]

        // Laminar
        Nu1 = 3.66 + 1.2*a_r^(-0.8); // Nu for thermally developed flow
        Nu11 =  fg*(Re*Pr*Dh/Ltube)^(1/3); // Nu for thermally developing flow
        Nu111 = (2/(1 + 22*Pr))^(1/6)*(Re*Pr*Dh/Ltube)^(1/2);

        // Transition
        Nu11_Lam_Tran =  fg*(2300*Pr*Dh/Ltube)^(1/3); // Nu for thermally developing flow
        Nu111_Lam_Tran = (2/(1 + 22*Pr))^(1/6)*(2300*Pr*Dh/Ltube)^(1/2);
        Nu_Lam_Tran = (Nu1^3 + Nu11_Lam_Tran^3 + Nu111_Lam_Tran^3)^(1/3); // Mean Nu for hydrodinamic developing and thermally developing flow
        k_Tur_Tran = 1.07 + 900/10e3 - 0.63/(1 + 10*Pr); // factor for turb. regime
        Re1_Tur_Tran = 10e3* ((1+a_r^2)*log(a_r) + (1 - a_r^2))/((1-a_r)^2*log(a_r));
        f_Tur_Tran = (1.8* log10(Re1_Tur_Tran) - 1.5)^(-2); // friction factor for turbulent regime
        Nu_Tur_Tran = f_corr*((f_Tur_Tran/8)*10e3*Pr)/(k_Tur_Tran + 12.7*(f_Tur_Tran/8)^(0.5)
                                                                                            *(Pr^(2/3) - 1))*(1+(Dh/Ltube)^(2/3));
        g = (Re - 2300)/(1e4 - 2300);

        // Turbulent
        k = 1.07 + 900/Re - 0.63/(1 + 10*Pr); // factor for turbulent regime
        Re_1 = Re* ((1+a_r^2)*log(a_r) + (1 - a_r^2))/((1-a_r)^2*log(a_r));
        f = (1.8 * log10(Re_1) - 1.5)^(-2); // friction factor (turbulent regime)
        f_corr = 0.75*a_r^(-0.17);     // Annular Tube Nu correction Factor

        if Re<2300 then
        Nu = (Nu1^3 + Nu11^3 + Nu111^3)^(1/3); // Mean Nu for hydrodinamic developing and thermally developing flow
        else if (Re>=2300 and Re<1e4) then
        Nu = (1-g)*Nu_Lam_Tran + g*Nu_Tur_Tran;

        else
        Nu = f_corr*((f/8)*Re*Pr)/(k + 12.7*(f/8)^(0.5)
                                                      *(Pr^(2/3) - 1))*(1+(Dh/Ltube)^(2/3)); // Modified Gnielinski, Ref. [5]
        end if;
          end if;

        h = k_c*Nu/Dh; // heat convection coefficient

        end Tube_In_Tube_to_validate_with_VDIAtalas;

        model Heat_exchangers

        // model that selects wich heat exchanger you want to use

          import SI = Modelica.SIunits;
          import Modelica.Constants.pi;

        type ShellSideMassVelocity = Real(quantity= "ShellSideMassVelocity", unit="kg/(s.m2)");

        /******************************  GRAPHICS  ************************************ */

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  BOOLEAN PARAMETERS  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        parameter Boolean NumberOfTubes_Given = false;
        outer parameter Boolean Bell_Delaware_Method = true;
        outer parameter Boolean Calculated_Volume
            "if true then the volume is calculated from L_input and d_input";
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  HEAT TRANSFER PARAMETERS  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        outer parameter SI.MassFlowRate mDot_c = 30 "Coolant mass flow rate";   // OUTER
        outer parameter Real Ft = 0.875
            "temperature factor (1 = countercurrent)";                             // minimum value is 0.75

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /******************************************  VARIABLES AND PARAMETERS FOR HEX1 AND HEX2   **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /******************************************  SHELL AND TUBE PARAMETERS AND VARIABLES  ********************************************************************************************************************************************************************************************************/
        /************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        // PARAMETERS (Inner in "Tank System")
        outer parameter SI.Length d_tank
            "tank diameter - used if Calculated_Volume=true";
        outer parameter SI.Volume V_input
            "Inner volume of cylindrical tank - used if Calculated_Volume=false";

        SI.Length Ds = (if Calculated_Volume == false then (4*V_input/(Modelica.Constants.pi*Ltube))^(1/2) else d_tank)
            "tank diameter";

        outer parameter SI.Length ID "inner tube diameter = 2*delta";
        outer parameter SI.Length Ltube = 1 "Tube length";
        outer parameter Real PR = 1.25
            "Tube pitch to tube diameter ratio: Pt/OD";                             // Tube pitch to tube outer diameter ratio (Pt/OD) http://www.engineeringpage.com/technology/thermal/pitch.html
        // outer parameter SI.Length Ds = 0.4 "Shell Diameter";
        outer parameter Modelica.SIunits.Conversions.NonSIunits.Angle_deg theta_tp = 30
            "Layout angle (30,60,45,90deg)";                                                                        // in deg
        outer parameter Integer Nt_input = 100
            "Used if NumberOfTubes_Given = true";
        outer parameter Real Nss = 1
            "Number of Sealing strips (pairs) in one baffle spacing";
        outer parameter Real B0=0.5
            "Initial value for the baffle space 0<B<1 - percentage (>0.2 (>5cm in general))";                       // Baffle space as percentage of the shell diameter "Ds"  // OUTER
        outer parameter Real Bc = 25e-2
            "Baffle cut (20-49%) with 20-25% optimum";
        outer parameter Integer Nb = 4
            "Number of baffles - input if NumberOfTubes_Given = true";
        outer parameter Real Lo_star = 1
            "1=<Lo_star<2 --> Ratio between the inlet baffle spacing and the central baffle spacing"; // Ratio between the outlet baffle spacing and the central baffle spacing
        outer parameter Real Li_star = 1
            "1=<Li_star<2 --> Ratio between the inlet baffle spacing and the central baffle spacing"; // Ratio between the inlet baffle spacing and the central baffle spacing
        outer parameter Integer N_shell = 1; // Number of shells

        // FINAL PARAMETERS
        final parameter SI.Length OD = ID/0.8 "inner tube diameter";
        final parameter SI.Length Pt = PR*OD; // tube pitch
        final parameter SI.Length Cl_shell=0.025;   /* [m] - Distance between shell and tube bundle. From: "Mechanical Constraints of thermal design of Shell and Tube Exchangers" Available at: https://www.academia.edu */
        final parameter SI.Length CTP=0.93;  // tube count constant for one tube pass [1]
        final parameter SI.Length Cl_tube = Pt - OD; // Distance between tubes
        final parameter SI.Length Ratio = Cl_tube./Pt;
        final parameter SI.Length Lpl=0;  //  =0 for single tube pass. Expresses the effect of the tube lane partition bypass width;
        final parameter SI.Length Ltb = 0.4e-3; // Diametral clearance between tube outside diameter and baffle hole

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  Constants and Constant exponents  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        constant Real a1_1 = 1.4;
        constant Real a1_2 = 1.36;
        constant Real a1_3 = 0.593;
        constant Real a1_4 = 0.321;
        constant Real a2_1 = -0.667;
        constant Real a2_2 = -0.657;
        constant Real a2_3 = -0.477;
        constant Real a2_4 = -0.388;
        constant Real g_c = 1; // [6]

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /*******************************************************  VARIABLES **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        // Geometry
        SI.Length Dn;  // Nozzles diameter from [6] Ch. 5 (See table 5.3)
        Real x;
        Real r_s;
        Real r_lm;
        Real f_s; // Friction factor [1]
        SI.Length Lts; // Tube sheet thickness
        SI.Length Lto; // Shell length
        SI.Length B_space; // Actual baffle space (recalculated)
        SI.Length Ltp_eff;
        SI.Area A_shell_BD;
        SI.Angle theta_ds; // in rad to be used in the formulae
        SI.Angle theta_ctl; // in rad to be used in the formulae
        Real Fw;  // Fraction of tubes in the baffle window
        Real Fc;  // Fraction of tubes in pure cross flow
        SI.Area Swg; // Gross window flow area
        SI.Area Swt; // Segmental baffle window area
        SI.Area Sw; // Net crossflow area throughout one baffle window
        SI.Length Dw;  // Equivalent hydraulic diameter of a segmental baffle window (used only for DeltaP calculations)
        SI.Length Lpp; // Effective tube row distance in the flow direction (function of theta_tp)
        Integer Nb_integer;  // number of baffles (rounded number = integer)
        Real Re_s; // Reynolds number (shell side)
        Real CL;   /* tube layout constant (value valid for inclination angle by 30deg and 60deg) =1 for 45 or 90deg => more tube density for 30 and 60deg. See Ref. [1] */
        Real B;     // Baffle space as percentage of the shell diameter "Ds"   (B=B_space/Ds)
        Real a;
        Real a1;
        Real a2;
        Real a3;
        Real a4;
        Real b1;
        Real b2;
        Real b3;
        Real b4;
        Real b;
        SI.Velocity vf_Kern; // for water and similar fluids it should be in the range: 0.6-1.5 m/s "Heat Exchanger Design Book"
        SI.Length De; // equivalent diameter
        SI.CoefficientOfHeatTransfer h_BD;
        SI.CoefficientOfHeatTransfer h_kern;
        SI.CoefficientOfHeatTransfer h_id;
        SI.CoefficientOfHeatTransfer h_shell_and_tube = (if Bell_Delaware_Method == false then h_kern else h_BD);

        SI.Area A_shell_kern;
        Integer Nt_max_shell_and_tube; // maximum number of tubes that fit in the shell
        Integer Nt_w; // Number of tubes in window area
        Real Nt_cc; // Number of tube rows in crossflow (i.e. between tbaffle tips)
        Real Nt_cw; // Effective number of tube rows crossed in the baffle window
        Integer Nc; //  total number of tube rows crossed in the entire heat exchanger
        SI.Area Ahex;
        SI.Area Sb; // Bypass area between the shell and the tube bundle within one baffle
        SI.Area Stb; // Tube to baffle hole leakagearea for one baffle (to calculate baffle leakage effect parameters J1 and R1)
        SI.Area Ssb; // Shell-to-baffle leakage area for one baffle (to calculate baffle leakage effect parameters J1 and R1)
        Real Fsbp; // Ratio of the bypass area (Sb) to the overall crossflow area (Sm) - used for J1 and R1
        ShellSideMassVelocity Gs_kern; // Shell side mass velocity
        ShellSideMassVelocity Gs_BD; // Shell side mass velocity
        ShellSideMassVelocity Gw; // Window mass velocity
        ShellSideMassVelocity Gn; // Nozzles mass velocity
        Real Nu_kern;
        Real r_ss; // =Nss/Nt_cc
        Real C_bh;
        Real C_bp;
        Real n_s; // exponent in "Js" (0.6 ia a value valid for turbulent flow)
        Real n_b; // exponent in "Js" (0.6 ia a value valid for turbulent flow)
        SI.Length B_space0; // [m] - Baffle space in absolute value to be used initially. the actual one is recalculate = B_space
        SI.Length Lbb;  // [m] Bundle to shell clearance
        SI.Length D_otl; // Tube bank OTL diameter
        SI.Length Lsb; // [7]
        SI.Length D_ctl; // Bundle diameter
        Real Nb_real;   // number of baffles (real number = decimal)

        // Properties
        Real Re_n; // Reynolds number in the nozzles
        Real Pr; // Prandtl number

        // Correction factors (J -> heat transfer; R -> Pressure losses)
        Real Ji; // Colburn factor
        Real J1; // Correction factor for baffle leakage effects, including both shell-to-baffle and tube-to-baffle leakage.
        Real Jc; // Correction factor for segmental baffle window (considers baffle cut and baffle spacing)
        Real Jr; // Only applicable for Laminar flow (Re_s<=100 and fullt effective for Re_s<20) with the maximum limit Jr=0.4
        Real Jb; // Correction factors for bundle bypass effects for heat transfer
        Real Js; // Heat transfer correction factor for unequal baffle spacing at inlet and/or outlet
        Real J_tot; // Cumulative correction factor (Good design: J_tot>=0.6)

        Real Rb; // Correction factors for bundle bypass effects for pressure drop
        Real R1;
        Real Rs; // Pressure drop correction factor for unequal baffle spacing at inlet and/or outlet

        // Pressure Losses
        SI.Pressure DeltaP_id; // Ideal shell-side tube-bank pressure drop
        SI.Pressure DeltaP_tot_no_nozzles; // Total pressure dorp on the shell side neglecting the nozzles
        SI.Pressure DeltaP_tot_with_nozzles; // Total pressure dorp on the shell side including the nozzles
        SI.Pressure DeltaP_c; // Pressure drop in all central baffle spaces
        SI.Pressure DeltaP_w; // Pressure drop in one window
        SI.Pressure DeltaP_e; // Pressure drop in the entrance and exit baffle spaces
        SI.Pressure DeltaP_w_tot; // Pressure drop in all windowst baffle spaces
        SI.Pressure DeltaP_n; // Pressure drop in the shell's nozzles

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /******************************************  TUBE IN TUBE PARAMETERS AND VARIABLES  ********************************************************************************************************************************************************************************************************/
        /************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        /************************************************  PARAMETERS  **********************************************************************************************************************************************************************************************************************************/

        outer parameter Real a_r "Aspect ratio; 0<a<1 (ID/OD)"; // outer
        final parameter SI.Length OD_tube_in_tube = ID/a_r;  // outer tube's inner diameter
        final parameter SI.Length Dh = OD_tube_in_tube-ID;  // hydraulic diameter

        final parameter SI.Length ODo = OD_tube_in_tube/0.8; // Outer tube's outer diameter (for Nt_max calculation) => considers tube's thickness
        final parameter Real fg = 1.615*(1 + 0.14*a_r^(-0.5));

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  VARIABLES  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        SI.MassFlowRate mDot_c_tube;  // mass flow rate per tube (in the annulus)

        final parameter SI.KinematicViscosity nu_c = mu_c/rho_c;

        Integer Nt_max_tube_in_tube;

        SI.Area A_cross_tube_in_tube; // Cross flow area
        Real Re_tube_in_tube;
        SI.Velocity v_c_tube_in_tube;

        // VARIABLES FOR Nu CALCULATION IN LAMINAR REGIME
        Real Nu1;  // Nu for thermally developed flow
        Real Nu11;  // Nu for thermally developing flow
        Real Nu111;

        // VARIABLES FOR Nu CALCULATION IN TRANSITION REGIME
        Real Nu11_Lam_Tran; // Nu for thermally developing flow
        Real Nu111_Lam_Tran;
        Real Nu_Lam_Tran; // Mean Nu for hydrodinamic developing and thermally developing flow

        Real k_Tur_Tran;  // factor for turb. regime
        Real Re1_Tur_Tran;
        Real f_Tur_Tran; // friction factor for turbulent regime
        Real Nu_Tur_Tran;
        Real g;

        //  VARIABLES FOR Nu CALCULATION IN TURBULENT REGIME
        Real k;  // factor for turbulent regime
        Real Re_1;
        Real f;  // friction factor (turbulent regime)
        Real f_corr;     // Annular Tube Nu correction Factor

        Real Nu;
        SI.CoefficientOfHeatTransfer h_tube_in_tube;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  PROPERTIES  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        outer parameter SI.Density rho_c;
        outer parameter SI.SpecificHeatCapacity cp_c;
        outer parameter SI.ThermalConductivity k_c;
        outer parameter SI.DynamicViscosity mu_c;

        constant SI.Density rho_w = 1000;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  HEAT EXCHANGER TYPE  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        // Heat exchanger type stored in records:
        inner replaceable parameter
            MHSS.Heat_Transfer.Heat_Exchangers.Heat_Exchanger_Configuration HEXCHANGER
            constrainedby
            MHSS.Heat_Transfer.Heat_Exchangers.Heat_Exchanger_Configuration                                     annotation (choicesAllMatching=true);
         final parameter Integer HEX = HEXCHANGER.HEX;  // HEX = 1 => shell_and_tube; HEX=2 => tube_in_tube

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  EQUATIONS  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

          Ports.HeatTransferCoefficientPort heatTransferCoefficientPort annotation (Placement(transformation(extent={{-6,32},
                    {14,52}}),                                                                                                    iconTransformation(extent={{-14,26},
                    {14,54}})));

        equation
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************************  SHELL AND TUBE CALCULATIONS ********************************************************************************************************************************************************************************************************/
        /************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        if NumberOfTubes_Given == false then
          Nt_max_shell_and_tube = integer(floor(0.785 * (CTP/CL)* Ds^2/(PR^2*OD^2) + 0.5)); // Eq. 8.9, Ref. [1]
          Nb_integer = if (Nb_real > 0) then integer(floor(Nb_real + 0.5)) else integer(ceil(Nb_real - 0.5));
        else
         Nt_max_shell_and_tube = Nt_input;
            Nb_integer = Nb;

        end if;

         // heatTransferCoefficientPort.h = h_shell_and_tube;
         // heatTransferCoefficientPort.Nt = Nt_max_shell_and_tube;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  KERN METHOD  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        A_shell_kern = (B0*Ds^2*Cl_tube)/Pt; /* [m2] Cross flow area shell-side; 
  modified equation (original Eq 8.15 pp. 309 [1]: 
  Ds*Cl_tube*Distance_Baffle/Pt, here rearranged with this substitution: 
  Distance_Baffle = B*Ds where: 20%< B <100% ). Note that Cl_tube/Pt =Ratio
  that is: Cl_tube/Pt = (Pt-OD)/Pt = 1-1/PR and thus is a constant.
  At the end for each simulation where Ds and B are constants, A_shell_kern is a
  constant as well. */

        vf_Kern = mDot_c/(rho_c*A_shell_kern);
        Ahex = Modelica.Constants.pi*Ltube*OD*Nt_max_shell_and_tube;
        Pr = (cp_c*mu_c)/k_c;
        Gs_kern = mDot_c/A_shell_kern;
        De = 4*(Pt^2 - Modelica.Constants.pi*OD^2/4)/(Modelica.Constants.pi*OD); // shell hydraulic (or equivalent) diameter - square layout
        Nu_kern = 0.36*(De*Gs_kern/mu_c)^0.55 * (cp_c*mu_c/k_c)^(1/3);
        h_kern = k_c/De*Nu_kern;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  BELL-DELAWARE METHOD **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        /**************************************************  GEOMETRY CALCULATIONS  ********************************************************************************************************************************************************************************************************************************************/

        B_space0= B0*Ds;
        Lbb = (12*1e-3 + 0.005*Ds); // [m] Bundle to shell clearance
        D_otl = Ds-Lbb; // Tube bank OTL diameter
        Lsb = (3.1*1e-3 + 0.004*Ds); // [7]
        Nb_real= Ltube/B_space0 - 1;  // number of baffles (real number = decimal)
        D_ctl = D_otl - OD; // Bundle diameter

        A_shell_BD = B_space*(Lbb + D_ctl/Ltp_eff*(Pt-OD));
        Lto = Ltube + 2*Lts;
        B_space = Ltube/(Nb_integer + 1);  // exact baffle spacing
        B = B_space/Ds;
        h_id = Ji*cp_c*Gs_BD*(1/Pr)^(2/3); //  the parameter that considers the different dynamic viscosity at the wall and in the bulk fluid according to the temperature is not considered here yet (compute it when you have the temperature)
        Re_s = OD*mDot_c/(mu_c*A_shell_BD); // used for the Bell-Delaware equation, Ref. [1]
        Gs_BD = mDot_c/A_shell_BD;

        theta_ds = 2*(acos(1-2*Bc)); // Here "B" is given already in decimal form, tehrefore do not write B/100 as in the original formula of the book
        theta_ctl = 2*(acos(Ds/D_ctl*(1-2*Bc)));

        Fw = theta_ctl/(2*pi) - sin(theta_ctl)/(2*pi);
        Fc = 1-2*Fw;

        Swg = pi/4*Ds^2*(theta_ds/(2*pi) - sin(theta_ds)/(2*pi));
        Swt = Nt_max_shell_and_tube*Fw*pi/4*OD^2;
        Sw = Swg-Swt;
        Sb = B_space*(Ds - D_otl + Lpl);
        Ssb = pi/2*Ds*Lsb*((2*pi-theta_ds)/(2*pi));
        Stb = pi/4*((OD + Ltb)^2 - OD^2)*Nt_max_shell_and_tube*(1-Fw); // [7, 8]

        Nt_w = integer(floor(Nt_max_shell_and_tube*Fw + 0.5));
        // Nt_cc = integer(floor(Ds/Lpp*(1-2*Bc) + 0.5));
        Nt_cc = Ds/Lpp*(1-2*Bc);
        // Nt_cw = integer(floor(0.8/Lpp*(Ds*Bc - (Ds - D_ctl)/2) + 0.5));
        Nt_cw = 0.8/Lpp*(Ds*Bc - (Ds - D_ctl)/2);
        Nc = integer(floor((Nt_cc+Nt_cw)*(Nb_integer+1) + 0.5));

        Fsbp= Sb/A_shell_BD;
        Dw = 4*Sw/(pi*OD*Nt_w + pi*Ds*theta_ds/(2*pi));

        if theta_tp == 30 then
        Lpp = 0.866*Pt;
        else if  (theta_tp == 60) then
         Lpp = 0.866*Pt;
             else if (theta_tp == 90) then
           Lpp = Pt;
               else
        Lpp = 0.707*Pt; //  value for theta_tp = 45deg
        end if;
        end if;
        end if;

        if Ds*0.1>25e-3 then
        Lts = 0.1*Ds;
        else
          Lts = 25.4e-3;
        end if;

        if (theta_tp == 30 or theta_tp ==90) then
        Ltp_eff = Pt;
        else if  (theta_tp == 60) then
         Ltp_eff = 0.866*Pt;
          else
        Ltp_eff = 0.707*Pt; //  value for theta_tp = 45deg
        end if;
        end if;

        x = (-0.15*(1+r_s) + 0.8);

        r_s = Ssb/(Ssb + Stb);

        r_lm = (Ssb + Stb)/A_shell_BD;

        r_ss = Nss/Nt_cc;

        f_s = b1*(1.33/(Pt/OD))^b*(Re_s)^(b2); // [1]

        Gn = mDot_c/(pi*Dn^2/4);

        Re_n = 4*mDot_c/(pi*Dn*mu_c);

        Gw = mDot_c/(A_shell_BD*Sw)^(1/2);

        /**************************************************  COEFFICIENTS AND EXPONENTS  ********************************************************************************************************************************************************************************************************************************************/

        if (theta_tp == 30 or theta_tp ==60) then
          a3 = 1.450;
          a4 = 0.519;
          b3 = 7.00;
          b4 = 0.5;
         CL = 0.87;
        if Re_s<10 then
        a1 = a1_1;
        a2 = a2_1;
        b1= 48;
        b2 = -1;
        else if Re_s >= 10 and Re_s < 100 then
        a1= a1_2;
        a2= a2_2;
        b1= 45.1;
        b2 = -0.973;

        else if Re_s>=100 and Re_s<1000 then
        a1 = a1_3;
        a2 = a2_3;
        b1= 4.570;
        b2 = -0.476;

             else if Re_s>=1000 and Re_s<10000 then
        a1 = a1_4;
        a2 = a2_4;
        b1= 0.486;
        b2 = -0.152;

             else
        a1 = a1_4;
        a2 = a2_4;
        b1= 0.372;
        b2 = -0.123;

        end if;
        end if;
        end if;
        end if;
        else if theta_tp == 45 then
          a3 = 1.930;
          a4 = 0.5;
          b3 = 6.59;
          b4 = 0.520;
        CL = 1;

            if Re_s<10 then
        a1 = 1.550;
        a2 = -0.667;
        b1= 32;
        b2 = -1;

        else if Re_s >= 10 and Re_s < 100 then
        a1= 0.498;
        a2= -0.656;
        b1= 26.2;
        b2 = -0.913;

        else if Re_s>=100 and Re_s<1000 then
        a1 = 0.730;
        a2 = -0.5;
        b1= 3.5;
        b2 = -0.476;

             else  if Re_s>=1000 and Re_s<10000 then
        a1 = 0.37;
        a2 = -0.396;
        b1= 0.333;
        b2 = -0.136;

             else
            a1 = 0.37;
        a2 = -0.396;
        b1= 0.303;
        b2 = -0.126;

        end if;
        end if;
        end if;
        end if;

        else    // so this is for the case theta_tp = 90deg
          a3 = 1.187;
          a4 = 0.370;
          b3 = 6.30;
          b4 = 0.378;
        CL = 1;

          if Re_s<10 then
        a1 = 0.970;
        a2 = -0.667;
        b1= 35;
        b2 = -1;

        else if Re_s >= 10 and Re_s < 100 then
        a1= 0.90;
        a2= -0.631;
        b1= 32.1;
        b2 = -0.963;

        else if Re_s>=100 and Re_s<1000 then
        a1 = 0.408;
        a2 = -0.46;
        b1= 6.09;
        b2 = -0.602;
             else if Re_s>=1000 and Re_s<10000 then
        a1 = 0.107;
        a2 = -0.266;
        b1= 0.0815;
        b2 = +0.022;

             else
        a1 = 0.370;
        a2 = -0.395;
        b1= 0.391;
        b2 = -0.148;

             end if;
        end if;
          end if;
          end if;
        end if;
        end if;

        a = a3/(1 + 0.14*(Re_s)^(a4));
        b = b3/(1 + 0.14*(Re_s)^b4);

        if r_ss< 0.5 then
        Rb = exp(-C_bp*Fsbp*(1-(2*r_ss)^(1/3)));
        else
        Rb = min(exp(-C_bp*Fsbp*(1-(2*r_ss)^(1/3))), 1);
        end if;

        if Re_s >= 100 then
        C_bh = 1.35;
        C_bp = 3.7;
        else
        C_bh = 1.25;
        C_bp = 4.5;
        end if;

        if Ds>=4*0.0254 and Ds<=10*0.0254 then
          Dn = 2*0.0254;
        elseif Ds>=12*0.0254 and Ds<=17.5*0.0254 then
          Dn = 3*0.0254;
          elseif Ds>=19.25*0.0254 and Ds<=21.25*0.0254 then
          Dn = 4*0.0254;
          elseif Ds>=23*0.0254 and Ds<=29*0.0254 then
          Dn = 6*0.0254;

          elseif Ds>=31*0.0254 and Ds<=37*0.0254 then
          Dn = 8*0.0254;

          elseif Ds>=39*0.0254 and Ds<=42*0.0254 then
          Dn = 10*0.0254;
        else
          Dn = 12*0.0254;

        end if;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  Heat Transfer and Pressure-drop correction factors  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        Ji = a1*(1.33/PR)^a*Re_s^(a2); // Colburn coefficient
        Jc = 0.55 + 0.72*Fc;
        J1 = 0.44*(1-r_s) + (1-0.44*(1-r_s))*Modelica.Constants.e^(-2.2*r_lm);
        Jb = min(exp(-C_bh*Fsbp*(1 - (2*r_ss)^(1/3))), 1);
        J_tot = Jc*J1*Jb*Js*Jr;

        if Re_s >= 100 then
        Jr = 1;
          n_s = 0.6;
          n_b = 0.2;
        Js = ((Nb_integer - 1) + Li_star^(1-n_s) + Lo_star^(1-n_s))/((Nb_integer - 1) + (Li_star) + (Lo_star)); // [6] by Serth
        DeltaP_w = (2+0.6*Nt_cw)*Gw^2/(2*g_c*rho_c);

        elseif (20<Re_s and Re_s<100) then
          Jr = max( 1.51/(Nc^(0.18)) + ((20-Re_s)/80)*(1.51/(Nc^(0.18)) - 1), 0.4); // maximum limit for Jr
           n_s = 1/3;
          n_b = 1;
        Js = (((Nb_integer - 1) + Li_star^(1-n_s) + Lo_star^(1-n_s))/((Nb_integer - 1) + (Li_star) + (Lo_star)) + 1) /2; // Foor laminar flow Js is about halfwas between 1 and the value calculated for turbulent conditions
        DeltaP_w = 26*Gw*mu_c/(g_c*rho_c)*(Nt_cw/(Pt-OD) + B_space/(Dw^2)) + 2*Gw^2/(g_c*rho_c);
        else
            Jr = max( 1.51/(Nc^(0.18)), 0.4); // maximum limit for Jr
             n_s = 1/3;
          n_b = 1;
        Js = (((Nb_integer - 1) + Li_star^(1-n_s) + Lo_star^(1-n_s))/((Nb_integer - 1) + (Li_star) + (Lo_star)) + 1) /2; // Foor laminar flow Js is about halfwas between 1 and the value calculated for turbulent conditions
        DeltaP_w = 26*Gw*mu_c/(g_c*rho_c)*(Nt_cw/(Pt-OD) + B_space/(Dw^2)) + 2*Gw^2/(g_c*rho_c);
        end if;

        R1 = Modelica.Constants.e^(-1.33*(1+r_s)*r_lm^(x));
        Rs =  ((1/Li_star)^(2-n_b) + (1/Lo_star)^(2-n_b)); // [7, 8]

        /***************************************************************  Heat Transfer coefficient (Bell_Delaware)  *******************************************************************************************************************************************************************************************************************************/
        h_BD = h_id*Jc*J1*Jb*Js*Jr;

         /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /***************************************************************  PRESSURE DROPS  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        DeltaP_id = (2*f_s*Nt_cc*(Gs_BD)^(2))/(g_c*rho_c); // Ideal tube bank pressure drop
        DeltaP_c = (Nb_integer -1)*(DeltaP_id*Rb*R1); // Combined pressure drop of all the interior crossflow sections
        DeltaP_e = 2*DeltaP_id*(1 + Nt_cw/Nt_cc)*Rb*Rs; // Combined pressure drop for the entrance and exit sections
        // DeltaP_w = pressure loss in windows depends on Re_s, thus is inserted in an "if-statement"
        DeltaP_w_tot = DeltaP_w*Nb_integer*R1; // Pressure drop in all windows
        DeltaP_tot_no_nozzles = DeltaP_e+ DeltaP_w_tot + DeltaP_c; // Total Pressure drop (neglecting pressure drop in the nozzles)

          if Re_n >= 100 then
            DeltaP_n = 2e-13*N_shell*Gn^2/(rho_c/rho_w);

        else
            DeltaP_n = 4e-13*N_shell*Gn^2/(rho_c/rho_w);
        end if;

        DeltaP_tot_with_nozzles = DeltaP_tot_no_nozzles + DeltaP_n;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /******************************************************************  TUBE IN TUBE CALCULATIONS  ********************************************************************************************************************************************************************************************************/
        /************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        A_cross_tube_in_tube = (pi/4*OD_tube_in_tube^2) - (pi/4*ID^2); // [m2] Cross flow area
        v_c_tube_in_tube = mDot_c_tube/(A_cross_tube_in_tube*rho_c); // coolant velocity
        Re_tube_in_tube = v_c_tube_in_tube*Dh/nu_c;           // Reynolds number

        mDot_c_tube = mDot_c/Nt_max_tube_in_tube; // mass flow rate per tube (in the annulus)
        Nt_max_tube_in_tube = integer(floor(0.785 * (CTP/CL)* Ds^2/(PR^2*ODo^2) + 0.5)); // Eq. 8.9, Ref. [1]

        // Laminar
        Nu1 = 3.66 + 1.2*a_r^(-0.8); // Nu for thermally developed flow
        Nu11 =  fg*(Re_tube_in_tube*Pr*Dh/Ltube)^(1/3); // Nu for thermally developing flow
        Nu111 = (2/(1 + 22*Pr))^(1/6)*(Re_tube_in_tube*Pr*Dh/Ltube)^(1/2);

        // Transition
        Nu11_Lam_Tran =  fg*(2300*Pr*Dh/Ltube)^(1/3); // Nu for thermally developing flow
        Nu111_Lam_Tran = (2/(1 + 22*Pr))^(1/6)*(2300*Pr*Dh/Ltube)^(1/2);
        Nu_Lam_Tran = (Nu1^3 + Nu11_Lam_Tran^3 + Nu111_Lam_Tran^3)^(1/3); // Mean Nu for hydrodinamic developing and thermally developing flow
        k_Tur_Tran = 1.07 + 900/10e3 - 0.63/(1 + 10*Pr); // factor for turb. regime
        Re1_Tur_Tran = 10e3* ((1+a_r^2)*log(a_r) + (1 - a_r^2))/((1-a_r)^2*log(a_r));
        f_Tur_Tran = (1.8* log10(Re1_Tur_Tran) - 1.5)^(-2); // friction factor for turbulent regime
        Nu_Tur_Tran = f_corr*((f_Tur_Tran/8)*10e3*Pr)/(k_Tur_Tran + 12.7*(f_Tur_Tran/8)^(0.5)
                                                                                            *(Pr^(2/3) - 1))*(1+(Dh/Ltube)^(2/3));
        g = (Re_tube_in_tube - 2300)/(1e4 - 2300);

        // Turbulent
        k = 1.07 + 900/Re_tube_in_tube - 0.63/(1 + 10*Pr); // factor for turbulent regime
        Re_1 = Re_tube_in_tube* ((1+a_r^2)*log(a_r) + (1 - a_r^2))/((1-a_r)^2*log(a_r));
        f = (1.8 * log10(Re_1) - 1.5)^(-2); // friction factor (turbulent regime)
        f_corr = 0.75*a_r^(-0.17);     // Annular Tube Nu correction Factor

        if Re_tube_in_tube<2300 then
        Nu = (Nu1^3 + Nu11^3 + Nu111^3)^(1/3); // Mean Nu for hydrodinamic developing and thermally developing flow
        else if (Re_tube_in_tube>=2300 and Re_tube_in_tube<1e4) then
        Nu = (1-g)*Nu_Lam_Tran + g*Nu_Tur_Tran;

        else
        Nu = f_corr*((f/8)*Re_tube_in_tube*Pr)/(k + 12.7*(f/8)^(0.5)
                                                                   *(Pr^(2/3) - 1))*(1+(Dh/Ltube)^(2/3)); // Modified Gnielinski, Ref. [5]
        end if;
          end if;

        h_tube_in_tube = k_c*Nu/Dh; // heat convection coefficient
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************************  SELECTION OF HEAT EXCHANGER  ********************************************************************************************************************************************************************************************************/
        /************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        if HEX == 1 then

          heatTransferCoefficientPort.h = h_shell_and_tube;
            heatTransferCoefficientPort.Nt = Nt_max_shell_and_tube;

        else if HEX == 2 then
            heatTransferCoefficientPort.h = h_tube_in_tube;
              heatTransferCoefficientPort.Nt = Nt_max_tube_in_tube;

        else
            heatTransferCoefficientPort.h = h_shell_and_tube;
            heatTransferCoefficientPort.Nt = Nt_max_shell_and_tube;
        end if;
        end if;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************ REFERENCES ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /*
[1] S. Kakac, H. Liu. "Heat Exchangers Selection, rating and Thermal 
    Design". 2nd Edition, CRC Press. ISBN: 0-8493-0902-6.
  
[2] Serth R.W. Process Heat Transfer, Principles and Applications. Ed. 2007 ISBN: 978-0-12-373588-1

[3] Heat Exchanger Design Book

[4] Wolverine Tube Heat Transfer Data Book
    
*/

          annotation (Icon(graphics={Bitmap(extent={{-90,-62},{94,56}},
                    fileName=
                      "modelica://HySDeP/Graphics/Shell_And_Tube__And_Tube_In_Tube.png")}));
        end Heat_exchangers;

        model Heat_exchangers_for1D_MODEL

        // model that selects wich heat exchanger you want to use

          import SI = Modelica.SIunits;
          import Modelica.Constants.pi;

        type ShellSideMassVelocity = Real(quantity= "ShellSideMassVelocity", unit="kg/(s.m2)");

        /******************************  GRAPHICS  ************************************ */

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  BOOLEAN PARAMETERS  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        parameter Boolean NumberOfTubes_Given = false;
        outer parameter Boolean Bell_Delaware_Method;
         outer parameter Boolean Calculated_Volume
            "if true then the volume is calculated from L_input and d_input";

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  HEAT TRANSFER PARAMETERS  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        outer parameter SI.MassFlowRate mDot_c "Coolant mass flow rate";
        outer parameter Real Ft = 0.875
            "temperature factor (1 = countercurrent)";                             // minimum value is 0.75

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /******************************************  VARIABLES AND PARAMETERS FOR HEX1 AND HEX2   **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /******************************************  SHELL AND TUBE PARAMETERS AND VARIABLES  ********************************************************************************************************************************************************************************************************/
        /************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        // PARAMETERS (Inner in "Tank System")
         outer parameter SI.Length d_tank
            "tank diameter - used if Calculated_Volume=true";

         outer parameter SI.Volume V_input
            "Inner volume of cylindrical tank - used if Calculated_Volume=false";

        SI.Length Ds = (if Calculated_Volume == false then (4*V_input/(Modelica.Constants.pi*Ltube))^(1/2) else d_tank)
            "tank diameter";

         outer parameter SI.Length ID "inner tube diameter = 2*delta";
         outer parameter SI.Length Ltube "Tube length";
         outer parameter Real PR = 1.25
            "Tube pitch to tube diameter ratio: Pt/OD";                              // Tube pitch to tube outer diameter ratio (Pt/OD) http://www.engineeringpage.com/technology/thermal/pitch.html
        // outer parameter SI.Length Ds = 0.4 "Shell Diameter";
         outer parameter Modelica.SIunits.Conversions.NonSIunits.Angle_deg theta_tp = 30
            "Layout angle (30,60,45,90deg)";                                                                // in deg
         outer parameter Integer Nt_input = 100
            "Used if NumberOfTubes_Given = true";
         outer parameter Real Nss = 1
            "Number of Sealing strips (pairs) in one baffle spacing";
         outer parameter Real B0=0.5
            "Initial value for the baffle space 0<B<1 - percentage (>0.2 (>5cm in general))";                      // Baffle space as percentage of the shell diameter "Ds"  // OUTER
         outer parameter Real Bc = 25e-2
            "Baffle cut (20-49%) with 20-25% optimum";
         outer parameter Integer Nb = 4
            "Number of baffles - input if NumberOfTubes_Given = true";
         outer parameter Real Lo_star = 1
            "1=<Lo_star<2 --> Ratio between the inlet baffle spacing and the central baffle spacing"; // Ratio between the outlet baffle spacing and the central baffle spacing
         outer parameter Real Li_star = 1
            "1=<Li_star<2 --> Ratio between the inlet baffle spacing and the central baffle spacing"; // Ratio between the inlet baffle spacing and the central baffle spacing
         outer parameter Integer N_shell = 1; // Number of shells

        // FINAL PARAMETERS
        final parameter SI.Length OD = ID/0.8 "inner tube diameter";
        final parameter SI.Length Pt = PR*OD; // tube pitch
        final parameter SI.Length Cl_shell=0.025;   /* [m] - Distance between shell and tube bundle. From: "Mechanical Constraints of thermal design of Shell and Tube Exchangers" Available at: https://www.academia.edu */
        final parameter SI.Length CTP=0.93;  // tube count constant for one tube pass [1]
        final parameter SI.Length Cl_tube = Pt - OD; // Distance between tubes
        final parameter SI.Length Ratio = Cl_tube./Pt;
        final parameter SI.Length Lpl=0;  //  =0 for single tube pass. Expresses the effect of the tube lane partition bypass width;
        final parameter SI.Length Ltb = 0.4e-3; // Diametral clearance between tube outside diameter and baffle hole

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  Constants and Constant exponents  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        constant Real a1_1 = 1.4;
        constant Real a1_2 = 1.36;
        constant Real a1_3 = 0.593;
        constant Real a1_4 = 0.321;
        constant Real a2_1 = -0.667;
        constant Real a2_2 = -0.657;
        constant Real a2_3 = -0.477;
        constant Real a2_4 = -0.388;
        constant Real g_c = 1; // [6]

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /*******************************************************  VARIABLES **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        // Geometry
        SI.Length Dn;  // Nozzles diameter from [6] Ch. 5 (See table 5.3)
        Real x;
        Real r_s;
        Real r_lm;
        Real f_s; // Friction factor [1]
        SI.Length Lts; // Tube sheet thickness
        SI.Length Lto; // Shell length
        SI.Length B_space; // Actual baffle space (recalculated)
        SI.Length Ltp_eff;
        SI.Area A_shell_BD;
        SI.Angle theta_ds; // in rad to be used in the formulae
        SI.Angle theta_ctl; // in rad to be used in the formulae
        Real Fw;  // Fraction of tubes in the baffle window
        Real Fc;  // Fraction of tubes in pure cross flow
        SI.Area Swg; // Gross window flow area
        SI.Area Swt; // Segmental baffle window area
        SI.Area Sw; // Net crossflow area throughout one baffle window
        SI.Length Dw;  // Equivalent hydraulic diameter of a segmental baffle window (used only for DeltaP calculations)
        SI.Length Lpp; // Effective tube row distance in the flow direction (function of theta_tp)
        Integer Nb_integer;  // number of baffles (rounded number = integer)
        Real Re_s; // Reynolds number (shell side)
        Real CL;   /* tube layout constant (value valid for inclination angle by 30deg and 60deg) =1 for 45 or 90deg => more tube density for 30 and 60deg. See Ref. [1] */
        Real B;     // Baffle space as percentage of the shell diameter "Ds"   (B=B_space/Ds)
        Real a;
        Real a1;
        Real a2;
        Real a3;
        Real a4;
        Real b1;
        Real b2;
        Real b3;
        Real b4;
        Real b;
        SI.Velocity vf_Kern; // for water and similar fluids it should be in the range: 0.6-1.5 m/s "Heat Exchanger Design Book"
        SI.Length De; // equivalent diameter
        SI.CoefficientOfHeatTransfer h_BD;
        SI.CoefficientOfHeatTransfer h_kern;
        SI.CoefficientOfHeatTransfer h_id;
        SI.CoefficientOfHeatTransfer h_shell_and_tube = (if Bell_Delaware_Method == false then h_kern else h_BD);

        SI.Area A_shell_kern;
        Integer Nt_max_shell_and_tube; // maximum number of tubes that fit in the shell
        Integer Nt_w; // Number of tubes in window area
        Real Nt_cc; // Number of tube rows in crossflow (i.e. between tbaffle tips)
        Real Nt_cw; // Effective number of tube rows crossed in the baffle window
        Integer Nc; //  total number of tube rows crossed in the entire heat exchanger
        SI.Area Ahex;
        SI.Area Sb; // Bypass area between the shell and the tube bundle within one baffle
        SI.Area Stb; // Tube to baffle hole leakagearea for one baffle (to calculate baffle leakage effect parameters J1 and R1)
        SI.Area Ssb; // Shell-to-baffle leakage area for one baffle (to calculate baffle leakage effect parameters J1 and R1)
        Real Fsbp; // Ratio of the bypass area (Sb) to the overall crossflow area (Sm) - used for J1 and R1
        ShellSideMassVelocity Gs_kern; // Shell side mass velocity
        ShellSideMassVelocity Gs_BD; // Shell side mass velocity
        ShellSideMassVelocity Gw; // Window mass velocity
        ShellSideMassVelocity Gn; // Nozzles mass velocity
        Real Nu_kern;
        Real r_ss; // =Nss/Nt_cc
        Real C_bh;
        Real C_bp;
        Real n_s; // exponent in "Js" (0.6 ia a value valid for turbulent flow)
        Real n_b; // exponent in "Js" (0.6 ia a value valid for turbulent flow)
        SI.Length B_space0; // [m] - Baffle space in absolute value to be used initially. the actual one is recalculate = B_space
        SI.Length Lbb;  // [m] Bundle to shell clearance
        SI.Length D_otl; // Tube bank OTL diameter
        SI.Length Lsb; // [7]
        SI.Length D_ctl; // Bundle diameter
        Real Nb_real;   // number of baffles (real number = decimal)

        // Properties
        Real Re_n; // Reynolds number in the nozzles
        Real Pr; // Prandtl number

        // Correction factors (J -> heat transfer; R -> Pressure losses)
        Real Ji; // Colburn factor
        Real J1; // Correction factor for baffle leakage effects, including both shell-to-baffle and tube-to-baffle leakage.
        Real Jc; // Correction factor for segmental baffle window (considers baffle cut and baffle spacing)
        Real Jr; // Only applicable for Laminar flow (Re_s<=100 and fullt effective for Re_s<20) with the maximum limit Jr=0.4
        Real Jb; // Correction factors for bundle bypass effects for heat transfer
        Real Js; // Heat transfer correction factor for unequal baffle spacing at inlet and/or outlet
        Real J_tot; // Cumulative correction factor (Good design: J_tot>=0.6)

        Real Rb; // Correction factors for bundle bypass effects for pressure drop
        Real R1;
        Real Rs; // Pressure drop correction factor for unequal baffle spacing at inlet and/or outlet

        // Pressure Losses
        SI.Pressure DeltaP_id; // Ideal shell-side tube-bank pressure drop
        SI.Pressure DeltaP_tot_no_nozzles; // Total pressure dorp on the shell side neglecting the nozzles
        SI.Pressure DeltaP_tot_with_nozzles; // Total pressure dorp on the shell side including the nozzles
        SI.Pressure DeltaP_c; // Pressure drop in all central baffle spaces
        SI.Pressure DeltaP_w; // Pressure drop in one window
        SI.Pressure DeltaP_e; // Pressure drop in the entrance and exit baffle spaces
        SI.Pressure DeltaP_w_tot; // Pressure drop in all windowst baffle spaces
        SI.Pressure DeltaP_n; // Pressure drop in the shell's nozzles

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /******************************************  TUBE IN TUBE PARAMETERS AND VARIABLES  ********************************************************************************************************************************************************************************************************/
        /************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        /************************************************  PARAMETERS  **********************************************************************************************************************************************************************************************************************************/

         outer parameter Real a_r = 0.5 "Aspect ratio; 0<a<1 (ID/OD)";
        final parameter SI.Length OD_tube_in_tube = ID/a_r;  // outer tube's inner diameter
        final parameter SI.Length Dh = OD_tube_in_tube-ID;  // hydraulic diameter

        final parameter SI.Length ODo = OD_tube_in_tube/0.8; // Outer tube's outer diameter (for Nt_max calculation) => considers tube's thickness
        final parameter Real fg = 1.615*(1 + 0.14*a_r^(-0.5));

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  VARIABLES  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        SI.MassFlowRate mDot_c_tube;  // mass flow rate per tube (in the annulus)

        final parameter SI.KinematicViscosity nu_c = mu_c/rho_c;

        Integer Nt_max_tube_in_tube;

        SI.Area A_cross_tube_in_tube; // Cross flow area
        Real Re_tube_in_tube;
        SI.Velocity v_c_tube_in_tube;

        // VARIABLES FOR Nu CALCULATION IN LAMINAR REGIME
        Real Nu1;  // Nu for thermally developed flow
        Real Nu11;  // Nu for thermally developing flow
        Real Nu111;

        // VARIABLES FOR Nu CALCULATION IN TRANSITION REGIME
        Real Nu11_Lam_Tran; // Nu for thermally developing flow
        Real Nu111_Lam_Tran;
        Real Nu_Lam_Tran; // Mean Nu for hydrodinamic developing and thermally developing flow

        Real k_Tur_Tran;  // factor for turb. regime
        Real Re1_Tur_Tran;
        Real f_Tur_Tran; // friction factor for turbulent regime
        Real Nu_Tur_Tran;
        Real g;

        //  VARIABLES FOR Nu CALCULATION IN TURBULENT REGIME
        Real k;  // factor for turbulent regime
        Real Re_1;
        Real f;  // friction factor (turbulent regime)
        Real f_corr;     // Annular Tube Nu correction Factor

        Real Nu;
        SI.CoefficientOfHeatTransfer h_tube_in_tube;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  PROPERTIES  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        outer parameter SI.Density rho_c;  // OUTER
        outer parameter SI.SpecificHeatCapacity cp_c;  // OUTER
        outer parameter SI.ThermalConductivity k_c;  // OUTER
        outer parameter SI.DynamicViscosity mu_c;  // OUTER

        constant SI.Density rho_w = 1000;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  HEAT EXCHANGER TYPE  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        // Heat exchanger type stored in records:
        inner replaceable parameter
            MHSS.Heat_Transfer.Heat_Exchangers.Heat_Exchanger_Configuration HEXCHANGER
            constrainedby
            MHSS.Heat_Transfer.Heat_Exchangers.Heat_Exchanger_Configuration                                     annotation (choicesAllMatching=true);
         final parameter Integer HEX = HEXCHANGER.HEX;  // HEX = 1 => shell_and_tube; HEX=2 => tube_in_tube

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  EQUATIONS  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

          Ports.Nt_h Nt_h
            annotation (Placement(transformation(extent={{-20,44},{0,64}})));
          Ports.Nt_Port Nt_Port
            annotation (Placement(transformation(extent={{-88,48},{-68,68}}),
                iconTransformation(extent={{-88,44},{-68,64}})));
        equation
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************************  SHELL AND TUBE CALCULATIONS ********************************************************************************************************************************************************************************************************/
        /************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        if NumberOfTubes_Given == false then
          Nt_max_shell_and_tube = integer(floor(0.785 * (CTP/CL)* Ds^2/(PR^2*OD^2) + 0.5)); // Eq. 8.9, Ref. [1]
          Nb_integer = if (Nb_real > 0) then integer(floor(Nb_real + 0.5)) else integer(ceil(Nb_real - 0.5));
        else
         Nt_max_shell_and_tube = Nt_input;
            Nb_integer = Nb;

        end if;

         // heatTransferCoefficientPort.h = h_shell_and_tube;
         // heatTransferCoefficientPort.Nt = Nt_max_shell_and_tube;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  KERN METHOD  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        A_shell_kern = (B0*Ds^2*Cl_tube)/Pt; /* [m2] Cross flow area shell-side; 
  modified equation (original Eq 8.15 pp. 309 [1]: 
  Ds*Cl_tube*Distance_Baffle/Pt, here rearranged with this substitution: 
  Distance_Baffle = B*Ds where: 20%< B <100% ). Note that Cl_tube/Pt =Ratio
  that is: Cl_tube/Pt = (Pt-OD)/Pt = 1-1/PR and thus is a constant.
  At the end for each simulation where Ds and B are constants, A_shell_kern is a
  constant as well. */

        vf_Kern = mDot_c/(rho_c*A_shell_kern);
        Ahex = Modelica.Constants.pi*Ltube*OD*Nt_max_shell_and_tube;
        Pr = (cp_c*mu_c)/k_c;
        Gs_kern = mDot_c/A_shell_kern;
        De = 4*(Pt^2 - Modelica.Constants.pi*OD^2/4)/(Modelica.Constants.pi*OD); // shell hydraulic (or equivalent) diameter - square layout
        Nu_kern = 0.36*(De*Gs_kern/mu_c)^0.55 * (cp_c*mu_c/k_c)^(1/3);
        h_kern = k_c/De*Nu_kern;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  BELL-DELAWARE METHOD **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        /**************************************************  GEOMETRY CALCULATIONS  ********************************************************************************************************************************************************************************************************************************************/

        B_space0= B0*Ds;
        Lbb = (12*1e-3 + 0.005*Ds); // [m] Bundle to shell clearance
        D_otl = Ds-Lbb; // Tube bank OTL diameter
        Lsb = (3.1*1e-3 + 0.004*Ds); // [7]
        Nb_real= Ltube/B_space0 - 1;  // number of baffles (real number = decimal)
        D_ctl = D_otl - OD; // Bundle diameter

        A_shell_BD = B_space*(Lbb + D_ctl/Ltp_eff*(Pt-OD));
        Lto = Ltube + 2*Lts;
        B_space = Ltube/(Nb_integer + 1);  // exact baffle spacing
        B = B_space/Ds;
        h_id = Ji*cp_c*Gs_BD*(1/Pr)^(2/3); //  the parameter that considers the different dynamic viscosity at the wall and in the bulk fluid according to the temperature is not considered here yet (compute it when you have the temperature)
        Re_s = OD*mDot_c/(mu_c*A_shell_BD); // used for the Bell-Delaware equation, Ref. [1]
        Gs_BD = mDot_c/A_shell_BD;

        theta_ds = 2*(acos(1-2*Bc)); // Here "B" is given already in decimal form, tehrefore do not write B/100 as in the original formula of the book
        theta_ctl = 2*(acos(Ds/D_ctl*(1-2*Bc)));

        Fw = theta_ctl/(2*pi) - sin(theta_ctl)/(2*pi);
        Fc = 1-2*Fw;

        Swg = pi/4*Ds^2*(theta_ds/(2*pi) - sin(theta_ds)/(2*pi));
        Swt = Nt_max_shell_and_tube*Fw*pi/4*OD^2;
        Sw = Swg-Swt;
        Sb = B_space*(Ds - D_otl + Lpl);
        Ssb = pi/2*Ds*Lsb*((2*pi-theta_ds)/(2*pi));
        Stb = pi/4*((OD + Ltb)^2 - OD^2)*Nt_max_shell_and_tube*(1-Fw); // [7, 8]

        Nt_w = integer(floor(Nt_max_shell_and_tube*Fw + 0.5));
        // Nt_cc = integer(floor(Ds/Lpp*(1-2*Bc) + 0.5));
        Nt_cc = Ds/Lpp*(1-2*Bc);
        // Nt_cw = integer(floor(0.8/Lpp*(Ds*Bc - (Ds - D_ctl)/2) + 0.5));
        Nt_cw = 0.8/Lpp*(Ds*Bc - (Ds - D_ctl)/2);
        Nc = integer(floor((Nt_cc+Nt_cw)*(Nb_integer+1) + 0.5));

        Fsbp= Sb/A_shell_BD;
        Dw = 4*Sw/(pi*OD*Nt_w + pi*Ds*theta_ds/(2*pi));

        if theta_tp == 30 then
        Lpp = 0.866*Pt;
        else if  (theta_tp == 60) then
         Lpp = 0.866*Pt;
             else if (theta_tp == 90) then
           Lpp = Pt;
               else
        Lpp = 0.707*Pt; //  value for theta_tp = 45deg
        end if;
        end if;
        end if;

        if Ds*0.1>25e-3 then
        Lts = 0.1*Ds;
        else
          Lts = 25.4e-3;
        end if;

        if (theta_tp == 30 or theta_tp ==90) then
        Ltp_eff = Pt;
        else if  (theta_tp == 60) then
         Ltp_eff = 0.866*Pt;
          else
        Ltp_eff = 0.707*Pt; //  value for theta_tp = 45deg
        end if;
        end if;

        x = (-0.15*(1+r_s) + 0.8);

        r_s = Ssb/(Ssb + Stb);

        r_lm = (Ssb + Stb)/A_shell_BD;

        r_ss = Nss/Nt_cc;

        f_s = b1*(1.33/(Pt/OD))^b*(Re_s)^(b2); // [1]

        Gn = mDot_c/(pi*Dn^2/4);

        Re_n = 4*mDot_c/(pi*Dn*mu_c);

        Gw = mDot_c/(A_shell_BD*Sw)^(1/2);

        /**************************************************  COEFFICIENTS AND EXPONENTS  ********************************************************************************************************************************************************************************************************************************************/

        if (theta_tp == 30 or theta_tp ==60) then
          a3 = 1.450;
          a4 = 0.519;
          b3 = 7.00;
          b4 = 0.5;
         CL = 0.87;
        if Re_s<10 then
        a1 = a1_1;
        a2 = a2_1;
        b1= 48;
        b2 = -1;
        else if Re_s >= 10 and Re_s < 100 then
        a1= a1_2;
        a2= a2_2;
        b1= 45.1;
        b2 = -0.973;

        else if Re_s>=100 and Re_s<1000 then
        a1 = a1_3;
        a2 = a2_3;
        b1= 4.570;
        b2 = -0.476;

             else if Re_s>=1000 and Re_s<10000 then
        a1 = a1_4;
        a2 = a2_4;
        b1= 0.486;
        b2 = -0.152;

             else
        a1 = a1_4;
        a2 = a2_4;
        b1= 0.372;
        b2 = -0.123;

        end if;
        end if;
        end if;
        end if;
        else if theta_tp == 45 then
          a3 = 1.930;
          a4 = 0.5;
          b3 = 6.59;
          b4 = 0.520;
        CL = 1;

            if Re_s<10 then
        a1 = 1.550;
        a2 = -0.667;
        b1= 32;
        b2 = -1;

        else if Re_s >= 10 and Re_s < 100 then
        a1= 0.498;
        a2= -0.656;
        b1= 26.2;
        b2 = -0.913;

        else if Re_s>=100 and Re_s<1000 then
        a1 = 0.730;
        a2 = -0.5;
        b1= 3.5;
        b2 = -0.476;

             else  if Re_s>=1000 and Re_s<10000 then
        a1 = 0.37;
        a2 = -0.396;
        b1= 0.333;
        b2 = -0.136;

             else
            a1 = 0.37;
        a2 = -0.396;
        b1= 0.303;
        b2 = -0.126;

        end if;
        end if;
        end if;
        end if;

        else    // so this is for the case theta_tp = 90deg
          a3 = 1.187;
          a4 = 0.370;
          b3 = 6.30;
          b4 = 0.378;
        CL = 1;

          if Re_s<10 then
        a1 = 0.970;
        a2 = -0.667;
        b1= 35;
        b2 = -1;

        else if Re_s >= 10 and Re_s < 100 then
        a1= 0.90;
        a2= -0.631;
        b1= 32.1;
        b2 = -0.963;

        else if Re_s>=100 and Re_s<1000 then
        a1 = 0.408;
        a2 = -0.46;
        b1= 6.09;
        b2 = -0.602;
             else if Re_s>=1000 and Re_s<10000 then
        a1 = 0.107;
        a2 = -0.266;
        b1= 0.0815;
        b2 = +0.022;

             else
        a1 = 0.370;
        a2 = -0.395;
        b1= 0.391;
        b2 = -0.148;

             end if;
        end if;
          end if;
          end if;
        end if;
        end if;

        a = a3/(1 + 0.14*(Re_s)^(a4));
        b = b3/(1 + 0.14*(Re_s)^b4);

        if r_ss< 0.5 then
        Rb = exp(-C_bp*Fsbp*(1-(2*r_ss)^(1/3)));
        else
        Rb = min(exp(-C_bp*Fsbp*(1-(2*r_ss)^(1/3))), 1);
        end if;

        if Re_s >= 100 then
        C_bh = 1.35;
        C_bp = 3.7;
        else
        C_bh = 1.25;
        C_bp = 4.5;
        end if;

        if Ds>=4*0.0254 and Ds<=10*0.0254 then
          Dn = 2*0.0254;
        elseif Ds>=12*0.0254 and Ds<=17.5*0.0254 then
          Dn = 3*0.0254;
          elseif Ds>=19.25*0.0254 and Ds<=21.25*0.0254 then
          Dn = 4*0.0254;
          elseif Ds>=23*0.0254 and Ds<=29*0.0254 then
          Dn = 6*0.0254;

          elseif Ds>=31*0.0254 and Ds<=37*0.0254 then
          Dn = 8*0.0254;

          elseif Ds>=39*0.0254 and Ds<=42*0.0254 then
          Dn = 10*0.0254;
        else
          Dn = 12*0.0254;

        end if;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************  Heat Transfer and Pressure-drop correction factors  **********************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        Ji = a1*(1.33/PR)^a*Re_s^(a2); // Colburn coefficient
        Jc = 0.55 + 0.72*Fc;
        J1 = 0.44*(1-r_s) + (1-0.44*(1-r_s))*Modelica.Constants.e^(-2.2*r_lm);
        Jb = min(exp(-C_bh*Fsbp*(1 - (2*r_ss)^(1/3))), 1);
        J_tot = Jc*J1*Jb*Js*Jr;

        if Re_s >= 100 then
        Jr = 1;
          n_s = 0.6;
          n_b = 0.2;
        Js = ((Nb_integer - 1) + Li_star^(1-n_s) + Lo_star^(1-n_s))/((Nb_integer - 1) + (Li_star) + (Lo_star)); // [6] by Serth
        DeltaP_w = (2+0.6*Nt_cw)*Gw^2/(2*g_c*rho_c);

        elseif (20<Re_s and Re_s<100) then
          Jr = max( 1.51/(Nc^(0.18)) + ((20-Re_s)/80)*(1.51/(Nc^(0.18)) - 1), 0.4); // maximum limit for Jr
           n_s = 1/3;
          n_b = 1;
        Js = (((Nb_integer - 1) + Li_star^(1-n_s) + Lo_star^(1-n_s))/((Nb_integer - 1) + (Li_star) + (Lo_star)) + 1) /2; // Foor laminar flow Js is about halfwas between 1 and the value calculated for turbulent conditions
        DeltaP_w = 26*Gw*mu_c/(g_c*rho_c)*(Nt_cw/(Pt-OD) + B_space/(Dw^2)) + 2*Gw^2/(g_c*rho_c);
        else
            Jr = max( 1.51/(Nc^(0.18)), 0.4); // maximum limit for Jr
             n_s = 1/3;
          n_b = 1;
        Js = (((Nb_integer - 1) + Li_star^(1-n_s) + Lo_star^(1-n_s))/((Nb_integer - 1) + (Li_star) + (Lo_star)) + 1) /2; // Foor laminar flow Js is about halfwas between 1 and the value calculated for turbulent conditions
        DeltaP_w = 26*Gw*mu_c/(g_c*rho_c)*(Nt_cw/(Pt-OD) + B_space/(Dw^2)) + 2*Gw^2/(g_c*rho_c);
        end if;

        R1 = Modelica.Constants.e^(-1.33*(1+r_s)*r_lm^(x));
        Rs =  ((1/Li_star)^(2-n_b) + (1/Lo_star)^(2-n_b)); // [7, 8]

        /***************************************************************  Heat Transfer coefficient (Bell_Delaware)  *******************************************************************************************************************************************************************************************************************************/
        h_BD = h_id*Jc*J1*Jb*Js*Jr;

         /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /***************************************************************  PRESSURE DROPS  ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        DeltaP_id = (2*f_s*Nt_cc*(Gs_BD)^(2))/(g_c*rho_c); // Ideal tube bank pressure drop
        DeltaP_c = (Nb_integer -1)*(DeltaP_id*Rb*R1); // Combined pressure drop of all the interior crossflow sections
        DeltaP_e = 2*DeltaP_id*(1 + Nt_cw/Nt_cc)*Rb*Rs; // Combined pressure drop for the entrance and exit sections
        // DeltaP_w = pressure loss in windows depends on Re_s, thus is inserted in an "if-statement"
        DeltaP_w_tot = DeltaP_w*Nb_integer*R1; // Pressure drop in all windows
        DeltaP_tot_no_nozzles = DeltaP_e+ DeltaP_w_tot + DeltaP_c; // Total Pressure drop (neglecting pressure drop in the nozzles)

          if Re_n >= 100 then
            DeltaP_n = 2e-13*N_shell*Gn^2/(rho_c/rho_w);

        else
            DeltaP_n = 4e-13*N_shell*Gn^2/(rho_c/rho_w);
        end if;

        DeltaP_tot_with_nozzles = DeltaP_tot_no_nozzles + DeltaP_n;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /******************************************************************  TUBE IN TUBE CALCULATIONS  ********************************************************************************************************************************************************************************************************/
        /************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        A_cross_tube_in_tube = (pi/4*OD_tube_in_tube^2) - (pi/4*ID^2); // [m2] Cross flow area
        v_c_tube_in_tube = mDot_c_tube/(A_cross_tube_in_tube*rho_c); // coolant velocity
        Re_tube_in_tube = v_c_tube_in_tube*Dh/nu_c;           // Reynolds number

        mDot_c_tube = mDot_c/Nt_max_tube_in_tube; // mass flow rate per tube (in the annulus)
        Nt_max_tube_in_tube = integer(floor(0.785 * (CTP/CL)* Ds^2/(PR^2*ODo^2) + 0.5)); // Eq. 8.9, Ref. [1]

        // Laminar
        Nu1 = 3.66 + 1.2*a_r^(-0.8); // Nu for thermally developed flow
        Nu11 =  fg*(Re_tube_in_tube*Pr*Dh/Ltube)^(1/3); // Nu for thermally developing flow
        Nu111 = (2/(1 + 22*Pr))^(1/6)*(Re_tube_in_tube*Pr*Dh/Ltube)^(1/2);

        // Transition
        Nu11_Lam_Tran =  fg*(2300*Pr*Dh/Ltube)^(1/3); // Nu for thermally developing flow
        Nu111_Lam_Tran = (2/(1 + 22*Pr))^(1/6)*(2300*Pr*Dh/Ltube)^(1/2);
        Nu_Lam_Tran = (Nu1^3 + Nu11_Lam_Tran^3 + Nu111_Lam_Tran^3)^(1/3); // Mean Nu for hydrodinamic developing and thermally developing flow
        k_Tur_Tran = 1.07 + 900/10e3 - 0.63/(1 + 10*Pr); // factor for turb. regime
        Re1_Tur_Tran = 10e3* ((1+a_r^2)*log(a_r) + (1 - a_r^2))/((1-a_r)^2*log(a_r));
        f_Tur_Tran = (1.8* log10(Re1_Tur_Tran) - 1.5)^(-2); // friction factor for turbulent regime
        Nu_Tur_Tran = f_corr*((f_Tur_Tran/8)*10e3*Pr)/(k_Tur_Tran + 12.7*(f_Tur_Tran/8)^(0.5)
                                                                                            *(Pr^(2/3) - 1))*(1+(Dh/Ltube)^(2/3));
        g = (Re_tube_in_tube - 2300)/(1e4 - 2300);

        // Turbulent
        k = 1.07 + 900/Re_tube_in_tube - 0.63/(1 + 10*Pr); // factor for turbulent regime
        Re_1 = Re_tube_in_tube* ((1+a_r^2)*log(a_r) + (1 - a_r^2))/((1-a_r)^2*log(a_r));
        f = (1.8 * log10(Re_1) - 1.5)^(-2); // friction factor (turbulent regime)
        f_corr = 0.75*a_r^(-0.17);     // Annular Tube Nu correction Factor

        if Re_tube_in_tube<2300 then
        Nu = (Nu1^3 + Nu11^3 + Nu111^3)^(1/3); // Mean Nu for hydrodinamic developing and thermally developing flow
        else if (Re_tube_in_tube>=2300 and Re_tube_in_tube<1e4) then
        Nu = (1-g)*Nu_Lam_Tran + g*Nu_Tur_Tran;

        else
        Nu = f_corr*((f/8)*Re_tube_in_tube*Pr)/(k + 12.7*(f/8)^(0.5)
                                                                   *(Pr^(2/3) - 1))*(1+(Dh/Ltube)^(2/3)); // Modified Gnielinski, Ref. [5]
        end if;
          end if;

        h_tube_in_tube = k_c*Nu/Dh; // heat convection coefficient
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************************  SELECTION OF HEAT EXCHANGER  ********************************************************************************************************************************************************************************************************/
        /************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/

        if HEX == 1 then

          Nt_h.h = h_shell_and_tube;
            Nt_h.Nt = Nt_max_shell_and_tube;
        Nt_Port.Nt = Nt_max_shell_and_tube;
        else if HEX == 2 then
            Nt_h.h = h_tube_in_tube;
              Nt_h.Nt = Nt_max_tube_in_tube;
              Nt_Port.Nt  = Nt_max_tube_in_tube;

        else
            Nt_h.h = h_shell_and_tube;
            Nt_h.Nt = Nt_max_shell_and_tube;
              Nt_Port.Nt = Nt_max_shell_and_tube;

        end if;
        end if;

        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /************************************************ REFERENCES ************************************************************************************************************************************************************************************************************************************/
        /**********************************************************************************************************************************************************************************************************************************************************************************************/
        /*
[1] S. Kakac, H. Liu. "Heat Exchangers Selection, rating and Thermal 
    Design". 2nd Edition, CRC Press. ISBN: 0-8493-0902-6.
  
[2] Serth R.W. Process Heat Transfer, Principles and Applications. Ed. 2007 ISBN: 978-0-12-373588-1

[3] Heat Exchanger Design Book

[4] Wolverine Tube Heat Transfer Data Book
    
*/

          annotation (Icon(graphics={Bitmap(extent={{-90,-62},{94,56}},
                    fileName=
                      "modelica://HySDeP/Graphics/Shell_And_Tube__And_Tube_In_Tube1.png")}),
              Diagram(graphics));
        end Heat_exchangers_for1D_MODEL;
      end Heat_Exchangers;
    end Heat_Transfer;

    package Kinetics

    end Kinetics;

    package TESTS

      model TEST_0D_MHSS

        Tanks_With_HEX.Zero_Dimensional_Models.Tank_system_with_Hexs
          tank_system_with_Hexs(
          VariableProperties=false,
          Calculated_Volume=true,
          a_r=0.5,
          Charging=true,
          mDot_c=20,
          redeclare Properties.MH.Ti1_1CrMn MH_properties,
          redeclare Properties.Coolant.Dex_Cool Coolant,
          redeclare Properties.Tube.Aluminum6061 Tube_properties,
          redeclare Heat_Transfer.Heat_Exchangers.SHELL_AND_TUBE
            HEXCHANGER_USED,
          PR=1.25,
          delta=0.1e-2,
          Tc_in=273.15)
          annotation (Placement(transformation(extent={{-72,-68},{16,26}})));

        Sources.Pressure_ramp_and_Tamb pressure_ramp_and_Tamb(
          duration=60,
          p0=101325,
          height=30000000)
          annotation (Placement(transformation(extent={{-56,-24},{76,124}})));
      equation
        connect(pressure_ramp_and_Tamb.p_T, tank_system_with_Hexs.p_T)
          annotation (Line(
            points={{11.32,8.56},{11.32,-23.14},{11.6,-23.14},{11.6,-21}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(pressure_ramp_and_Tamb.p_ramp, tank_system_with_Hexs.p_ramp)
          annotation (Line(
            points={{-57.32,91.44},{-57.32,21.98},{-62.32,21.98},{-62.32,-21}},
            color={0,0,0},
            smooth=Smooth.None));
        annotation (Diagram(graphics));
      end TEST_0D_MHSS;

      model TEST_1D_MHSS

        MHSS.Tanks_With_HEX.One_Dimensional_Models.Tank_system_1D tank_system_1D_New(
          eps=0.6,
          VariableProperties=true,
          Charging=true,
          mDot_c=20,
          redeclare HySDeP.MHSS.Properties.MH.Ti1_1CrMn MH_properties,
          redeclare HySDeP.MHSS.Properties.Coolant.Dex_Cool Coolant,
          redeclare HySDeP.MHSS.Properties.Tube.Aluminum6061 Tube_properties,
          a_r=0.5,
          redeclare Heat_Transfer.Heat_Exchangers.SHELL_AND_TUBE
            HEXCHANGER_USED,
          N=10,
          delta=2e-2,
          Tc_in=273.15)
          annotation (Placement(transformation(extent={{-80,-102},{36,-2}})));

        Sources.Pressure_ramp_and_Tamb pressure_ramp_and_Tamb(
          duration=60,
          p0=100000,
          height=30000000)
          annotation (Placement(transformation(extent={{-82,-34},{76,124}})));
      equation
        connect(pressure_ramp_and_Tamb.p_T, tank_system_1D_New.p_T) annotation (
            Line(
            points={{-3,0.76},{-3,-4.41},{26.72,-4.41},{26.72,-52}},
            color={0,0,0},
            smooth=Smooth.None));
        connect(pressure_ramp_and_Tamb.p_ramp, tank_system_1D_New.p_ramp)
          annotation (Line(
            points={{-82,89.24},{-82,0.01},{-71.88,0.01},{-71.88,-52}},
            color={0,0,0},
            smooth=Smooth.None));
        annotation (Diagram(graphics));
      end TEST_1D_MHSS;

      model TEST_1D_MHSS_1

        MHSS.Tanks_With_HEX.One_Dimensional_Models.Tank_system_1D tank_system_1D_New(
          eps=0.6,
          VariableProperties=true,
          Charging=true,
          mDot_c=20,
          redeclare HySDeP.MHSS.Properties.MH.Ti1_1CrMn MH_properties,
          redeclare HySDeP.MHSS.Properties.Coolant.Dex_Cool Coolant,
          redeclare HySDeP.MHSS.Properties.Tube.Aluminum6061 Tube_properties,
          a_r=0.5,
          redeclare Heat_Transfer.Heat_Exchangers.SHELL_AND_TUBE
            HEXCHANGER_USED,
          N=10,
          delta=2e-2,
          Tc_in=273.15)
          annotation (Placement(transformation(extent={{-22,10},{6,38}})));

        Sources.Pressure_ramp_and_Tamb pressure_ramp_and_Tamb(
          duration=60,
          p0=100000,
          height=30000000)
          annotation (Placement(transformation(extent={{-26,18},{46,90}})));
      equation
        connect(pressure_ramp_and_Tamb.p_ramp, tank_system_1D_New.p_ramp)
          annotation (Line(points={{-26.72,74.16},{-26.72,24},{-20.04,24}},
              color={0,0,0}));
        connect(pressure_ramp_and_Tamb.p_T, tank_system_1D_New.p_T) annotation (
           Line(points={{10.72,33.84},{10.72,24},{3.76,24}}, color={0,0,0}));
      end TEST_1D_MHSS_1;
    end TESTS;
    annotation (Icon(graphics={Bitmap(extent={{-78,-78},{82,66}}, fileName=
                "modelica://HySDeP/Graphics/MH_Tank.png")}), Diagram(graphics={
            Bitmap(extent={{-74,-66},{70,80}}, fileName=
                "modelica://HySDeP/Graphics/MH_Tank.png")}));
  end MHSS;

  package CHG
    model SAE_J2601 "Sets the parameters for the fueling based on J2601"
    //Decideds the average pressure ramp rate, final pressure in HSS,
    //final state of charge and temperature out of the fueling station based on J2601
    //using ambient temperature, starting pressure in HSS and refueling protocol as input
      import SI = Modelica.SIunits;
      /*************************************************/
      /*General parameters*/
      /*************************************************/
    parameter Integer Fueling_protocol=1 "Fueling protocol" annotation (choices(
    choice=1 "A70 1-7 kg",
     choice=2 "A70 7-10 kg",
    choice=3 "A35",
    choice=4 "B70 1-7 kg",
     choice=5 "B70 7-10 kg",
     choice=6 "B35",
     choice=7 "C35",
     choice=8 "D35"));
    parameter SI.Temperature T_amb=293 "Ambient temperature";
    parameter SI.Pressure P_amb=101000 "Ambient pressure";
    parameter SI.Pressure P_start=2e6 "Start pressure i the HSS";
      /*************************************************/
      /*APRR*/
      /*************************************************/

      Real T=T_amb;
     inner Real P=P_start;
      Real  T_cool;
    inner Real APRR;

      SI.Pressure FP;
      Real SOC;
     inner SI.Pressure P_ref;
      Real APRR_1;
      SI.Pressure FP_1;
      Real SOC_1;
      Real APRR_2;
      SI.Pressure FP_2;
      Real SOC_2;
      Real APRR_3;
      SI.Pressure FP_3;
      Real SOC_3;
      Real APRR_4;
      SI.Pressure FP_4;
      Real SOC_4;
      Real APRR_5;
      SI.Pressure FP_5;
      Real SOC_5;
      Real APRR_6;
      SI.Pressure FP_6;
      Real SOC_6;
      Real APRR_7;
      SI.Pressure FP_7;
      Real SOC_7;
      Real APRR_8;
      SI.Pressure FP_8;
      Real SOC_8;

      /*************************************************/
      /*Consstants from SAE*/
      /*************************************************/
    constant SI.Temperature T1=273.15-40;
    constant SI.Temperature T2=273.15-20;
    constant SI.Temperature T3=273.15;

    constant SI.Pressure P_ref1=70e6;
    constant SI.Pressure P_ref2=35e6;

      /*************************************************/
      /*Look up tables for SAE*/
      /*************************************************/
      // A70 1-7 kg
      // For the APRR (Average pressure ramprate)
      Modelica.Blocks.Tables.CombiTable1Ds APRR1(
        tableOnFile=true,
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
        tableName="APRR_A1",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/APRR.mat")
        annotation (Placement(transformation(extent={{-10,30},{10,50}})));

      // For FP (final pressure)
      Modelica.Blocks.Tables.CombiTable2D FP1(
        tableOnFile=true,
        tableName="FP_A1",
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
        fileName=
            "C:/Users/edro/Documents/Dymola/External files/Lookuptables/FP.mat")
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      //For SOC (state of charge)
      Modelica.Blocks.Tables.CombiTable2D SOC1(
        tableOnFile=true,
        tableName="SOC_A1",
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
        fileName=
            "C:/Users/edro/Documents/Dymola/External files/Lookuptables/SOC.mat")
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

      // A70 7-10 kg
      // For the APRR (Average pressure ramprate)
      Modelica.Blocks.Tables.CombiTable1Ds APRR2(
        tableOnFile=true,
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
        tableName="APRR_A2",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/APRR.mat");

      // For FP (final pressure)
      Modelica.Blocks.Tables.CombiTable2D FP2(
        tableOnFile=true,
        tableName="FP_A2",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/FP.mat",
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);

      //For SOC (state of charge)
      Modelica.Blocks.Tables.CombiTable2D SOC2(
        tableOnFile=true,
        tableName="SOC_A2",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/SOC.mat",
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);

      // A35
      // For the APRR (Average pressure ramprate)
      Modelica.Blocks.Tables.CombiTable1Ds APRR3(
        tableOnFile=true,
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
        tableName="APRR_A3",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/APRR.mat");

      // For FP (final pressure)
      Modelica.Blocks.Tables.CombiTable2D FP3(
        tableOnFile=true,
        tableName="FP_A3",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/FP.mat",
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);

      //For SOC (state of charge)
      Modelica.Blocks.Tables.CombiTable2D SOC3(
        tableOnFile=true,
        tableName="SOC_A3",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/SOC.mat",
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);

      // B70 1-7kg
      // For the APRR (Average pressure ramprate)
      Modelica.Blocks.Tables.CombiTable1Ds APRR4(
        tableOnFile=true,
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
        tableName="APRR_B1",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/APRR.mat");

      // For FP (final pressure)
      Modelica.Blocks.Tables.CombiTable2D FP4(
        tableOnFile=true,
        tableName="FP_B1",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/FP.mat",
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);

      //For SOC (state of charge)
      Modelica.Blocks.Tables.CombiTable2D SOC4(
        tableOnFile=true,
        tableName="SOC_B1",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/SOC.mat",
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);

      // B70 7-10kg
      // For the APRR (Average pressure ramprate)
      Modelica.Blocks.Tables.CombiTable1Ds APRR5(
        tableOnFile=true,
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
        tableName="APRR_B2",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/APRR.mat");

      // For FP (final pressure)
      Modelica.Blocks.Tables.CombiTable2D FP5(
        tableOnFile=true,
        tableName="FP_B2",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/FP.mat",
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);

      //For SOC (state of charge)
      Modelica.Blocks.Tables.CombiTable2D SOC5(
        tableOnFile=true,
        tableName="SOC_B2",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/SOC.mat",
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);

      // B35
      // For the APRR (Average pressure ramprate)
      Modelica.Blocks.Tables.CombiTable1Ds APRR6(
        tableOnFile=true,
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
        tableName="APRR_B3",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/APRR.mat");

      // For FP (final pressure)
      Modelica.Blocks.Tables.CombiTable2D FP6(
        tableOnFile=true,
        tableName="FP_B3",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/FP.mat",
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);

      //For SOC (state of charge)
      Modelica.Blocks.Tables.CombiTable2D SOC6(
        tableOnFile=true,
        tableName="SOC_B3",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/SOC.mat",
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);

      // C35
      // For the APRR (Average pressure ramprate)
      Modelica.Blocks.Tables.CombiTable1Ds APRR7(
        tableOnFile=true,
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
        tableName="APRR_C1",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/APRR.mat");

      // For FP (final pressure)
      Modelica.Blocks.Tables.CombiTable2D FP7(
        tableOnFile=true,
        tableName="FP_C1",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/FP.mat",
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);

      //For SOC (state of charge)
      Modelica.Blocks.Tables.CombiTable2D SOC7(
        tableOnFile=true,
        tableName="SOC_C1",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/SOC.mat",
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);

      // D35
      // For the APRR (Average pressure ramprate)
      Modelica.Blocks.Tables.CombiTable1Ds APRR8(
        tableOnFile=true,
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
        tableName="APRR_D1",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/APRR.mat");

      // For FP (final pressure)
      Modelica.Blocks.Tables.CombiTable2D FP8(
        tableOnFile=true,
        tableName="FP_D1",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/FP.mat",
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);

      //For SOC (state of charge)
      Modelica.Blocks.Tables.CombiTable2D SOC8(
        tableOnFile=true,
        tableName="SOC_D1",
        fileName="C:/Users/edro/Documents/Dymola/External files/Lookuptables/SOC.mat",
        smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);

      /*************************************************/
      /*Equations*/
      /*************************************************/
    equation
      FP1.u1 = T;
      FP1.u2 = P;
      FP1.y = FP_1;
      SOC1.u1 = T;
      SOC1.u2 = P;
      SOC1.y = SOC_1;
      APRR1.u = T;
      APRR1.y[1] = APRR_1;

      FP2.u1 = T;
      FP2.u2 = P;
      FP2.y = FP_2;
      SOC2.u1 = T;
      SOC2.u2 = P;
      SOC2.y = SOC_2;
      APRR2.u = T;
      APRR2.y[1] = APRR_2;

      FP3.u1 = T;
      FP3.u2 = P;
      FP3.y = FP_3;
      SOC3.u1 = T;
      SOC3.u2 = P;
      SOC3.y = SOC_3;
      APRR3.u = T;
      APRR3.y[1] = APRR_3;

      FP4.u1 = T;
      FP4.u2 = P;
      FP4.y = FP_4;
      SOC4.u1 = T;
      SOC4.u2 = P;
      SOC4.y = SOC_4;
      APRR4.u = T;
      APRR4.y[1] = APRR_4;

      FP5.u1 = T;
      FP5.u2 = P;
      FP5.y = FP_5;
      SOC5.u1 = T;
      SOC5.u2 = P;
      SOC5.y = SOC_5;
      APRR5.u = T;
      APRR5.y[1] = APRR_5;

      FP6.u1 = T;
      FP6.u2 = P;
      FP6.y = FP_6;
      SOC6.u1 = T;
      SOC6.u2 = P;
      SOC6.y = SOC_6;
      APRR6.u = T;
      APRR6.y[1] = APRR_6;

      FP7.u1 = T;
      FP7.u2 = P;
      FP7.y = FP_7;
      SOC7.u1 = T;
      SOC7.u2 = P;
      SOC7.y = SOC_7;
      APRR7.u = T;
      APRR7.y[1] = APRR_7;

      FP8.u1 = T;
      FP8.u2 = P;
      FP8.y = FP_8;
      SOC8.u1 = T;
      SOC8.u2 = P;
      SOC8.y = SOC_8;
      APRR8.u = T;
      APRR8.y[1] = APRR_8;

      /*************************************************/
      /*Deciding output values*/
      /*************************************************/

      if Fueling_protocol == 1 then //A70 1-7 kg
        APRR = APRR_1;
        FP = FP_1;
        SOC = SOC_1;
        T_cool=T1;
        P_ref=P_ref1;
      elseif Fueling_protocol == 2 then //A70 7-10 kg
        APRR = APRR_2;
        FP = FP_2;
        SOC= SOC_2;
        T_cool=T1;
        P_ref=P_ref1;
      elseif Fueling_protocol == 3 then //A35
        APRR = APRR_3;
        FP = FP_3;
        SOC = SOC_3;
        T_cool=T1;
        P_ref=P_ref2;
      elseif Fueling_protocol == 4 then //B70 1-7 kg
        APRR = APRR_4;
        FP = FP_4;
        SOC = SOC_4;
        T_cool=T2;
        P_ref=P_ref1;
      elseif Fueling_protocol == 5 then //B70 7-10 kg
        APRR = APRR_5;
        FP = FP_5;
        SOC = SOC_5;
        T_cool=T2;
        P_ref=P_ref1;
      elseif Fueling_protocol == 6 then //B35
        APRR = APRR_6;
        FP = FP_6;
        SOC = SOC_6;
        T_cool=T2;
        P_ref=P_ref2;
      elseif Fueling_protocol == 7 then //C35
        APRR = APRR_7;
        FP = FP_7;
        SOC = SOC_7;
        T_cool=T3;
        P_ref=P_ref2;
      elseif Fueling_protocol == 8 then //D35
        APRR = APRR_8;
        FP = FP_8;
        SOC = SOC_8;
        T_cool=T;
        P_ref=P_ref2;
      else
        APRR = 0;
        FP = 0;
        SOC = 0;
        T_cool=0;
        P_ref=0;
      end if;

      /*************************************************/
      /*Checking input values*/
      /*************************************************/

    algorithm
      assert(Fueling_protocol < 9, "Not a valid fueling procedure. Choose one between 1 and 8;
  1: A70 1-7 kg,
  2: A70 7-10 kg,
  3: A35,
  4: B70 1-7 kg,
  5: B70 7-10 kg,
  6: B35,
  7: C35,
  8: D35");
      assert(if Fueling_protocol == 3 or Fueling_protocol == 6 or
      Fueling_protocol == 7 or Fueling_protocol == 8 then P <= 35e6
         else P <= 70e6,
        "The initial HSS pressure is above allowance of fuelling procedure -
     No need for refilling, HSS is already full!");

      annotation (preferedView="text", Icon(graphics={Bitmap(
              extent={{-70,-62},{74,58}}, fileName=
                  "modelica://HySDeP/Graphics/SAE_J2601.png")}));
    end SAE_J2601;

    package Tanks "The different tanks used for hydrogen refueling"
      model Tank "Tank with editable Volume"
        import SI = Modelica.SIunits;

       /****************** Thermodynamic property call *************************/
      /* replaceable package Medium = CoolProp2Modelica.Media.Hydrogen (onePhase=true)
   constrainedby Modelica.Media.Interfaces.PartialMedium
                                               annotation (choicesAllMatching=true);
*/

       replaceable package Medium =
          ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
          Modelica.Media.Interfaces.PartialMedium;

             Medium.ThermodynamicState medium;
             Medium.ThermodynamicState mediumStream;
          //   Medium.ThermodynamicState medium1;
            // Medium.ThermodynamicState medium2;

       /******************** Connectors *****************************/

        Ports.FlowPort portA "port connection component to other components"
          annotation (Placement(transformation(extent={{50,-10},{70,10}},    rotation=
                 0), iconTransformation(extent={{18,-8},{36,8}})));
      //( m_flow(final start=m_flowStart))
        Ports.HeatFlow2 heatFlow "connection to heat transfer model"
          annotation (Placement(transformation(extent={{-112,-62},{-82,-32}}),
              iconTransformation(extent={{-10,-48},{10,-30}})));
        Ports.PressurePort pp "Connection for control of system"
          annotation (Placement(transformation(extent={{-114,38},{-72,78}}),
              iconTransformation(extent={{-10,30},{10,50}})));

       /****************** General parameters *******************/
        parameter Boolean Adiabatic = true "If true, adiabatic tank model" annotation(Dialog(group="Tank data"));
       outer SI.Volume  V "inner tank volume";

       /******************  Initial and start values *******************/

        //parameter SI.Pressure pInitial=1.013e5 "Initial pressure in the tank"
        parameter Boolean fixedInitialPressure = true "Fixed intial pressure"
                                  annotation(Dialog(group="Initial Values"));

        parameter SI.Temperature TInitial=T_amb
          "Initial temperature in the tank"
          annotation(Dialog(group="Initial Values"));

        parameter SI.MassFlowRate m_flowStart=0 "Initial mass flow rate"
          annotation(Dialog(group="Initial Values"));

        outer parameter SI.Temperature T_amb "Ambient temperature";

       /****************** variables *******************/

        SI.Mass M "Gas mass in control volume";
        Real drhodt;
        SI.Heat Q;
        SI.InternalEnergy U;
        SI.SpecificInternalEnergy u;
        SI.Pressure p;
        SI.SpecificEnthalpy h;
        Real HeatOfCompression "heat of compression";
        constant Real Counter=0;

      //Exergy varialbles
      //outer SI.Pressure P_amb;
      outer SI.SpecificEntropy s_0;
      outer SI.SpecificEnthalpy h_0;
      outer SI.SpecificInternalEnergy u_0;
      SI.Heat E_tank;
      SI.Heat E_stream;
      SI.SpecificEnthalpy e_tank;
      SI.SpecificEnthalpy e_stream;

      SI.Heat E_D;

      /****************** Initial equations *******************/
      initial equation
      h=Medium.specificEnthalpy_pT(p, TInitial);

      /****************** equations *******************/
      equation

      medium=Medium.setState_ph(p, h);
      u=(h-p*1/medium.d);
      U=u*M;
      HeatOfCompression=V*der(p);

      if Adiabatic == false then
      der(Q)=heatFlow.Q;
      else
      der(Q)=0;
      end if;

       der(h) = 1/M*(noEvent(actualStream(portA.h_outflow))*portA.m_flow - portA.m_flow*h + V*der(p)+der(Q))
          "Energy balance";

        p = portA.p;
        portA.h_outflow = h;
        M = V*medium.d "Mass in cv";
       drhodt = Medium.density_derp_h(medium)*der(p)+Medium.density_derh_p(medium)*der(h)
          "Derivative of density";

       drhodt*V = portA.m_flow "Mass balance";

      // Exergy
      if portA.m_flow >= 0 then
      mediumStream=Medium.setState_ph(p, inStream(portA.h_outflow));
      else
       mediumStream=Medium.setState_ph(p, h);
      end if;

      //u_0=h_0-P_amb*V;
      e_stream=mediumStream.h-h_0-T_amb*(mediumStream.s-s_0);
      e_tank=u-u_0-T_amb*(medium.s-s_0);
      der(E_tank)=e_tank*abs(portA.m_flow);
      der(E_stream)=e_stream*abs(portA.m_flow);
      E_tank=der(Q)*(1-T_amb/medium.T)+E_stream-E_D;
      //der(E_tank)=der(Q)*(1-T_amb/medium.T)+E_stream-E_D;

      //heatFlow.Q=der(Q);
      heatFlow.m_flow=portA.m_flow;
      heatFlow.T=medium.T;
      heatFlow.P=p;
      heatFlow.Counter=Counter;

      pp.p=p;

      /*
medium1 = Medium.setState_pT(700e5, );
medium2 = Medium.setState_pT(700e5, h);
*/
         annotation(Dialog(group="Tank data"),
                     Dialog(group="Initial Values"),
                    preferedView="text",Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-40},{120,40}},
              initialScale=0.1), graphics={Bitmap(
                extent={{-30,36},{30,-38}},
                imageSource=
                    "iVBORw0KGgoAAAANSUhEUgAABrIAAA9vCAIAAADZtbdkAAAAAXNSR0IArs4c6QAAAARnQU1BAACxjwv8YQUAAAAJcEhZcwAAXEYAAFxGARSUQ0EAAP+lSURBVHhe7N0/SOvrvuD/zWUyEwIDYSAQZEALmYmVXrBQGESQgRSDiDBos0BstLBwNSKMIDKgDBbaBIsUprlYWg0pLaawlFtZpkyZMuXvl3v8nH33vnvttX3W8k+eb16v6p57ztl7r4/J9/s87+P3+/zy/wEAUES9Xu8hxe3t7env3d3dxb/3G/1+P/4GAADkTBYEAPgg/X4/0trDw93dXbS309OdnZ3VV5iZmfllXI3+2eKf8je2t7fjT/gbNzc3MYLfeHp6ihkBAPBRZEEAgB8XWevhodvtRvc6PT0+Po4wtrq6sLAQ5YxE5XI5hvg3zWYz5nt6GkN/eIgfAwAA6WRBAIB/9fz8HMHp4eHq6ioq1Olps9mMOrW6Wi6XI1wxNpaWll5+OoeHhy8/svv7+5efo6eeAQC+SRYEACbFSyR6iX2/fW63Xq9HW6LQGo3Gy0989NN/SYe/PtH8/PwcnxIAgIkhCwIAhfL4+Pjw8HBzc3N6erq/v7/qMV5SVKvVl3S4sbHxkg4vLi5e0uHooxUfMgCAQpAFAYD8/LH9LS0tRdeB9zf6yL286/AlGjoyBQDIkSwIAIypp6enh4eH29vb09PTw8PDv/0K12pUGRg/L79p+HL+8svjyb1eLz7NAADjRxYEAD7TyxEf2h8FNjMzM/pU7+/vjz7kLwehDAaD+AIAAHweWRAAeHeDweDh4eHu7u709PT4+Pil/TnPlwn3cnry6Bsx+l6MviAj8YUBAPgQsiAA8MZ6vd7Dw8PFxcXLW//kP3i90fdl9K15eXHh1dXV6KvkxYUAwDuRBQGAn/L8/Nztdk9PT3d2dlY9/wvvo16vj75fo2/Z6Lt2e3v78PDQ7/fjSwgA8ENkQQAgwdPT0/39/enp6cbGhsN/4dM1Go3V1dXRV/Lu7u7x8TG+qAAAryALAgB/6uHvh4E0m81GoxEdAhhjCwsLL6ch+41CAOD7ZEEA4F8Mh8OHh4ebm5uXI0FmZmaiMQA5e3lZ4eHh4ejbPfqOxxceAEAWBIDJ1O/3Hx4erq6uDg8PV1dXq9VqJASg6GZmZl6ONLm/v39+fo6LAgAweWRBACi+l6OBT09PHQ0M/NHosrCzs3NxcTG6UAwGg7hwAABFJwsCQNG8PA58enq6vb092u3Hvh/gdV5OPT4+Pr69vXWMCQAUmCwIAEXw9PR0c3Ozs7OzsLAQO3uAN/LrMSbdbrfX68V1BwDInCwIAFka7czv7+9fjgeJjTvAh/j1GJOrq6uHh4fhcBgXJgAgK7IgAOTh10eDm81mvV6P3TnAGPjtMSb9fj8uWwDAeJMFAWB8PT4+vjwa3Gg0YvMNMPZmZma2t7dHl6+np6e4nAEA40cWBIAx0uv17u7uPBoMFEa5XH75RUKPGwPAuJEFAeAzeTQYmChLS0uHh4d3d3eeNQaATycLAsBH82gwwMjLs8ZXV1eeNQaATyELAsC7e3k0+PDw0KPBAN/0crqxZ40B4CPJggDw9gaDwa+PBler1dj1AvA6CwsLL88a93q9uLACAG9NFgSAt/H09HR1deXRYIC39euzxo+Pj3HBBQDegiwIAD+u1+vd3t6O9qtOCwH4AL8+a9ztdj1rDAA/SRYEgDSDweD+/n5/f99vBQJ8roWFhdHV2LPGAPBjZEEAeJWHh4fj4+OlpaXYjAIwTur1umeNASCJLAgAf+rldYHNZrNcLse+E4Cx9/Ks8cXFxegyHhd0AOAPZEEA+B2vCwQoktHFfGdnZ3Rh7/f7caEHAP5GFgQArwsEmAgLCwvHx8eeMgaAF7IgAJOr2+2O9oejXWLsFwGYDOVyeWNj4+bmxlklAEwyWRCAyfL09HRxcdFsNmNrCMBkazQah4eH9/f3w+EwbhUAMBlkQQCKr9fr3dzcbG9vV6vV2AUCwB80m82rq6vn5+e4fwBAocmCABTTYDC4u7vb39+fmZmJ3R4AvM7o3rGzszO6j4zuJnFfAYDCkQUBKBSvCwTgbS0tLV1cXDioBIDikQUByJ7XBQLwAarV6vb29u3tbb/fjzsQAORMFgQgSy/PCHtdIACfYmFh4fDwsNvtxm0JADIkCwKQk36/f3t7u7GxEdsyAPhU5XJ5dFe6ublxUAkA2ZEFAcjAy1HCq6ursQkDgPEzMzOzv79/f38/HA7jBgYAY0wWBGB8PT8/X1xcOD8EgOysrq6ObmFPT09xSwOA8SMLAjB2Rpuo4+PjRqMRWysAyFa9Xt/f3394eIibHACMDVkQgHEx2jIdHh7OzMzERgoACkQfBGDcyIIAfLL7+/vRNmm0WYptEwAUmj4IwJiQBQH4BMPh8O7ubnt7u1qtxiYJACaMPgjA55IFAfg4g8Hg9vZ2Y2OjXC7HlggAJl61Wt3Z2bm/v4/7JQB8CFkQgHfX7/dvbm6azWbsfgCAb9EHAfhIsiAA76XX611dXS0tLcVeB/hwU1NTc6+ztra2+RvNZjP+jd/z4D98DH0QgA8gCwLwxp6fn09PTxuNRuxsgN+r1WrR2H4f475+/XryOq1W65/GSafTiX+y39vb24s/2+/FH/73pqamYkDAb+iDALwfWRCAt/H09HR4eDgzMxP7GJgkUbbm5ubn5yN9bW5++fIl8tjJyeXlZSQ0Upydnb0McDTMl6kuLy+/jLpWq8X0YTLogwC8OVkQgJ/y8PCwv79fr9dj1wKF8Jpf6Ot0OtGu+FStVuvlJzL66bz8mNbX119+drOzs/EThQLRBwF4K7IgAMmGw+FoNzLak3jLGDkafW7nfvNrfUdHRy9RyS/0FdvLT3lka2vr5Ue/uLj4Ug9dysiUPgjAT5IFAUjQ7Xa3t7fL5XLsSGBcVSqVl+LzEoBefsvv+vo6EhF8y+gT8pIODw4OXj45vx69Eh8sGEsvffDu7m44HMYNGwBeQRYE4K89Pz8fHx97UphxUyqVXpLN+vr65ubmwcGB9sf7eTlZ5ejo6Ndc6JgUxk25XN7e3tYHAXglWRCAPzUYDG5vb5eWlmK3AZ/kpf01m81f258HfhkfLy83fDl2+eVEFE8l8+n0QQBeQxYE4Bu63e7Ozo6HhflIs7Ozv7a/vb29k5OT8/Pz6C6QodEHePQxfnmV4d/KtieR+QT6IADfIQsC8K96vd7x8fHMzExsJuAdTE1NvRz3sbu7e3JycnZ2FhEFJsDLk8gvhyZ7EpmPVC6Xd3Z2Hh4e4pYPALIgACPD4fD29nZ1dTW2DvB2pqenXyLgwcGBAgh/5uW0k98+iVypVOJbBG9qZmbm4uKi3+/HIgCACSYLAky0h4cHDwvzhqanpxcXFzf/dvKvR4Dh5/36JPL6+vqctxbypjY2Nu7u7mJBAMBEkgUBJlGv17u4uPCwMD9pdnZ2eXl5c3Pz6OjIGSDwMdrt9ksoXFtb875Cfl69Xj88PHx+fo4lAgCTRBYEmCAvDws3m83YCkCKubm5lZWVzc3NEwcBwzgZfR+Pjo5G383FxUVvKuSHLS0tjRYJg8EgFg0ATABZEGAiPD4+7u/ve/qMVyqVSnNzc2tray
              
              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"
                     +
                    "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AAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAEDYu3uYuBZsz9vdV1O6JTTSlEZCQkwAoyGACAKkwQkvkhMiC5HYE1hCTnDgwCSIwCPkBAdIYxLkgIBKZhw6JHRIiG5EyM0ICSvsd3d79bmn+9jneNl81Nr7eaLuvucD/2tX7dq/S9UGAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChO"
                     +
                    "FkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2Iz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",
                fileName=
                    "modelica://HydrogenRefuelingCoolPropAndrea/HydrogenTank.png"),
                Text(
                extent={{-24,14},{22,-14}},
                lineColor={0,0,0},
                textString="V")}), Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-40},{120,40}},
              initialScale=0.1), graphics),
              Documentation(info="<html>
        <a href=\"../Documentation/HydrogenLibaryDocumnetation.pdf\">PhD project by Erasmus Rothuizen</a><br><br>
        </html>"));
      end Tank;

      model Tank1 "Volume with one entrance"
        import SI = Modelica.SIunits;

       /****************** Thermodynamic property call *************************/
      /* replaceable package Medium = CoolProp2Modelica.Media.Hydrogen (onePhase=true)
   constrainedby Modelica.Media.Interfaces.PartialMedium
                                               annotation (choicesAllMatching=true);
*/

       replaceable package Medium =
          ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
          Modelica.Media.Interfaces.PartialMedium;

             Medium.ThermodynamicState medium;
             Medium.ThermodynamicState mediumStream;

       /******************** Connectors *****************************/

        Ports.FlowPort portA "port connection component to other components"
          annotation (Placement(transformation(extent={{50,-10},{70,10}},    rotation=
                 0), iconTransformation(extent={{20,-8},{38,8}})));
      //( m_flow(final start=m_flowStart))
        Ports.HeatFlow2 heatFlow "connection to heat transfer model"
          annotation (Placement(transformation(extent={{-112,-62},{-82,-32}}),
              iconTransformation(extent={{-10,-50},{10,-32}})));
        Ports.PressurePort pp "Connection for control of system"
          annotation (Placement(transformation(extent={{-114,38},{-72,78}}),
              iconTransformation(extent={{-10,30},{10,50}})));

       /****************** General parameters *******************/
        parameter Boolean Adiabatic = false "If true, adiabatic tank model" annotation(Dialog(group="Tank data"));
       parameter SI.Volume V=1 "Volume of the tank";

       /******************  Initial and start values *******************/

        //parameter SI.Pressure pInitial=1.013e5 "Initial pressure in the tank"
        parameter Boolean fixedInitialPressure = true "Fixed intial pressure"
                                  annotation(Dialog(group="Initial Values"));

        parameter SI.Temperature TInitial=T_amb
          "Initial temperature in the tank"
          annotation(Dialog(group="Initial Values"));

        parameter SI.MassFlowRate m_flowStart=0 "Initial mass flow rate"
          annotation(Dialog(group="Initial Values"));

        outer parameter SI.Temperature T_amb "Ambient temperature";

       /****************** variables *******************/

        SI.Mass M "Gas mass in control volume";
        Real drhodt;
        SI.Heat Q;
        SI.InternalEnergy U;
        SI.SpecificInternalEnergy u;
        SI.Pressure p;
        SI.SpecificEnthalpy h;
        Real HeatOfCompression "heat of compression";
        constant Real Counter=0;

      //Exergy varialbles
      //outer SI.Pressure P_amb;
      outer SI.SpecificEntropy s_0;
      outer SI.SpecificEnthalpy h_0;
      outer SI.SpecificInternalEnergy u_0;
      SI.Heat E_tank;
      SI.Heat E_stream;
      SI.SpecificEnthalpy e_tank;
      SI.SpecificEnthalpy e_stream;

      SI.Heat E_D;

      /****************** Initial equations *******************/
      initial equation
      h=Medium.specificEnthalpy_pT(p, TInitial);

      /****************** equations *******************/
      equation

      medium=Medium.setState_ph(p, h);
      u=(h-p*1/medium.d);
      U=u*M;
      HeatOfCompression=V*der(p);

      if Adiabatic == false then
      der(Q)=heatFlow.Q;
      else
      der(Q)=0;
      end if;

       der(h) = 1/M*(noEvent(actualStream(portA.h_outflow))*portA.m_flow - portA.m_flow*h
       + V*der(p)+der(Q)) "Energy balance";

        p = portA.p;
        portA.h_outflow = h;
        M = V*medium.d "Mass in cv";
       drhodt = Medium.density_derp_h(medium)*der(p)+Medium.density_derh_p(medium)*der(h)
          "Derivative of density";

       drhodt*V = portA.m_flow "Mass balance";

      // Exergy
      if portA.m_flow >= 0 then
      mediumStream=Medium.setState_ph(p, inStream(portA.h_outflow));
      else
       mediumStream=Medium.setState_ph(p, h);
      end if;

      //u_0=h_0-P_amb*V;
      e_stream=mediumStream.h-h_0-T_amb*(mediumStream.s-s_0);
      e_tank=u-u_0-T_amb*(medium.s-s_0);
      der(E_tank)=e_tank*abs(portA.m_flow);
      der(E_stream)=e_stream*abs(portA.m_flow);
      E_tank=der(Q)*(1-T_amb/medium.T)+E_stream-E_D;
      //der(E_tank)=der(Q)*(1-T_amb/medium.T)+E_stream-E_D;

      //heatFlow.Q=der(Q);
      heatFlow.m_flow=portA.m_flow;
      heatFlow.T=medium.T;
      heatFlow.P=p;
      heatFlow.Counter=Counter;

      pp.p=p;

         annotation(Dialog(group="Tank data"),
                     Dialog(group="Initial Values"),
                    preferedView="text",Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-40},{120,40}},
              initialScale=0.1), graphics={Bitmap(
                extent={{-36,34},{36,-34}},
                imageSource=
                    "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"
                     +
                    "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AAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAAAAAOdswUGbAQAAAMA5W3DQZgAAAABwzhYctBkAAAAAnLMFB20GAEDYu3uYuBZsz9vdV1O6JTTSlEZCQkwAoyGACAKkwQkvkhMiC5HYE1hCTnDgwCSIwCPkBAdIYxLkgIBKZhw6JHRIiG5EyM0ICSvsd3d79bmn+9jneNl81Nr7eaLuvucD/2tX7dq/S9UGAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChO"
                     +
                    "FkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2IzAAAAAChOFkyIzQAAAACgOFkwITYDAAAAgOJkwYTYDAAAAACKkwUTYjMAAAAAKE4WTIjNAAAAAKA4WTAhNgMAAACA4mTBhNgMAAAAAIqTBRNiMwAAAAAoThZMiM0AAAAAoDhZMCE2AwAAAIDiZMGE2AwAAAAAipMFE2Iz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",
                fileName=
                    "modelica://HydrogenRefuelingCoolPropAndrea/HydrogenTank.png")}),
                                   Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-40},{120,40}},
              initialScale=0.1), graphics),
              Documentation(info="<html>
        <a href=\"../Documentation/HydrogenLibaryDocumnetation.pdf\">PhD project by Erasmus Rothuizen</a><br><br>
        </html>"));
      end Tank1;

      model Tank2 "Volume with two entrances"
        import SI = Modelica.SIunits;
       /****************** Gas *************************/
          replaceable package Medium =
           CoolProp2Modelica.Interfaces.ExternalTwoPhaseMedium annotation(Dialog(group="Gas"));

             Medium.ThermodynamicState medium;
             Medium.ThermodynamicState mediumStreamA;
             Medium.ThermodynamicState mediumStreamB;

       /******************** Connectors *****************************/

        Ports.FlowPort portB(
                            m_flow(final start=m_flowStart))
          "port connection component to other components"
          annotation (Placement(transformation(extent={{50,-10},{70,10}},    rotation=
                 0), iconTransformation(extent={{92,-10},{114,10}})));

        Ports.HeatFlow2 heatFlow "connection to heat transfer model"
          annotation (Placement(transformation(extent={{-112,-62},{-82,-32}}),
              iconTransformation(extent={{-10,-48},{10,-30}})));
        Ports.PressurePort pp "Connection for control of system"
          annotation (Placement(transformation(extent={{-114,38},{-72,78}}),
              iconTransformation(extent={{-12,26},{12,50}})));
        Ports.FlowPort portA(m_flow(final start=m_flowStart))
          "port connection component to other components"
          annotation (Placement(transformation(extent={{50,-28},{70,-8}},    rotation=
                 0), iconTransformation(extent={{-104,-10},{-82,10}})));
       /****************** General parameters *******************/
        parameter Boolean Adiabatic = false "If true, adiabatic tank model" annotation(Dialog(group="Tank data"));
        parameter SI.Volume V=1 "Volume of the tank" annotation(Dialog(group="Tank data"));

       /******************  Initial and start values *******************/

        parameter SI.Pressure pInitial=1.013e5 "Initial pressure in the tank"
          annotation(Dialog(group="Initial Values"));
        parameter Boolean fixedInitialPressure = true "Fixed intial pressure"
                                  annotation(Dialog(group="Initial Values"));

        parameter SI.Temperature TInitial=T_amb
          "Initial temperature in the tank"
          annotation(Dialog(group="Initial Values"));

        parameter SI.MassFlowRate m_flowStart=0 "Initial mass flow rate"
          annotation(Dialog(group="Initial Values"));

       outer parameter SI.Temperature T_amb "Ambient temperature";

       /****************** variables *******************/

        SI.Mass M "Gas mass in control volume";
        Real drhodt;
        SI.HeatFlowRate Q;
        SI.InternalEnergy U;
        SI.SpecificInternalEnergy u;
        SI.Pressure p(start=pInitial, fixed=fixedInitialPressure);
        SI.SpecificEnthalpy h;
        Real HeatOfCompression "heat of compression";
        constant Real Counter=0;

      //Exergy varialbles
      //outer SI.Pressure P_amb;
      outer SI.SpecificEntropy s_0;
      // SI.SpecificEntropy s;
      outer SI.SpecificEnthalpy h_0;
      outer SI.SpecificInternalEnergy u_0;
      SI.Heat E_tank;
      SI.Heat E_streamA;
      SI.Heat E_streamB;
      SI.HeatFlowRate e_tank;
      SI.HeatFlowRate e_streamA;
      SI.HeatFlowRate e_streamB;
      SI.Heat E_D;

      /****************** Initial equations *******************/
      initial equation
      h=Medium.specificEnthalpy_pT(pInitial, TInitial);

      /****************** equations *******************/
      equation
      u=(h-p*1/medium.d);
      U=u*M;
      medium=Medium.setState_ph(p, h);

      HeatOfCompression=V*der(p);

      if Adiabatic == false then
      der(Q)=heatFlow.Q;
      else
      der(Q)=0;
      end if;

      der(h) = 1/M*(noEvent(actualStream(portA.h_outflow))*portA.m_flow - portA.m_flow*h
       + noEvent(actualStream(portB.h_outflow))*portB.m_flow - portB.m_flow*h
       + V*der(p)+der(Q));

        p = portA.p;
        portA.p-portB.p=0;

        portA.h_outflow = h;
        portB.h_outflow = h;

        M = V*medium.d "Mass in cv";
        drhodt = Medium.density_derp_h(medium)*der(p)+Medium.density_derh_p(medium)*der(h)
          "Derivative of density";

      drhodt*V = portA.m_flow + portB.m_flow "Mass balance";

      // Exergy
      if portA.m_flow >= 0 then
      mediumStreamA=Medium.setState_ph(p, inStream(portA.h_outflow));
      else
       mediumStreamA=Medium.setState_ph(p, h);
      end if;

      if portB.m_flow >= 0 then
      mediumStreamB=Medium.setState_ph(p, inStream(portA.h_outflow));
      else
       mediumStreamB=Medium.setState_ph(p, h);
      end if;

      //u_0=(h_0-P_amb*V);
      e_streamA=mediumStreamA.h-h_0-T_amb*(mediumStreamA.s-s_0);
      e_streamB=mediumStreamB.h-h_0-T_amb*(mediumStreamB.s-s_0);
      e_tank=u-u_0-T_amb*(medium.s-s_0);
      der(E_tank)=e_tank*(portA.m_flow+portB.m_flow) "Exergy in tank";
      der(E_streamA)=e_streamA*portA.m_flow "Exergy in stream from port A";
      der(E_streamB)=e_streamB*portB.m_flow "Exergy in stream from port B";
      E_tank=der(Q)*(1-T_amb/medium.T)+E_streamA+E_streamB-E_D "Exergy balance";

      //heatFlow.Q=der(Q);
      heatFlow.m_flow=portA.m_flow-portB.m_flow;
      heatFlow.T=medium.T;
      heatFlow.P=p;
      heatFlow.Counter=Counter;

      pp.p=p;

        annotation (preferedView="text",Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-40},{120,40}},
              initialScale=0.1), graphics={Bitmap(extent={{-96,42},{102,-42}},
                  fileName=
                    "modelica://HydrogenRefuelingCoolProp/Graphics/TankHorizontal.png")}),
                                   Diagram(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-40},{120,40}},
              initialScale=0.1), graphics),
              Documentation(info="<html>
        <a href=\"../Documentation/HydrogenLibaryDocumnetation.pdf\">PhD project by Erasmus Rothuizen</a><br><br>
        </html>"));
      end Tank2;

    end Tanks;

    package Ports "Flow connectors, 6 different flow connectors"

      connector FlowPort
        "Flow port passing on mass flow, pressure and enthalpy"
        import SI = Modelica.SIunits;
       SI.Pressure p "Pressure";
       flow SI.MassFlowRate m_flow "Mass flow rate";
       stream SI.SpecificEnthalpy h_outflow "Specific enthalpy";

        annotation (Icon(coordinateSystem(extent={{-80,-80},{80,80}}),
                         graphics={Ellipse(
                extent={{-80,80},{80,-80}},
                lineColor={0,0,0},
                lineThickness=0.5,
                fillColor={0,0,255},
                fillPattern=FillPattern.Solid)}),
                                  Diagram(coordinateSystem(extent={{-80,-80},{80,
                  80}}),                  graphics={Ellipse(
                extent={{-60,60},{60,-60}},
                lineColor={0,0,0},
                fillColor={0,0,255},
                fillPattern=FillPattern.Solid), Text(
                extent={{-100,94},{100,54}},
                lineColor={0,0,255},
                fillColor={255,153,0},
                fillPattern=FillPattern.Solid,
                textString=
                     "%name")}),
                     Documentation(info="<html>
        <a href=\"../Documentation/HydrogenLibaryDocumnetation.pdf\">PhD project by Erasmus Rothuizen</a><br><br>
        </html>"));
      end FlowPort;

      connector PressurePort "Passes on pressure"
        import SI = Modelica.SIunits;
        SI.Pressure p "Pressure";
        annotation (Icon(graphics={Ellipse(
                extent={{-80,80},{80,-80}},
                lineColor={0,0,0},
                lineThickness=0.5,
                fillColor={170,255,255},
                fillPattern=FillPattern.Solid)}),
                               Documentation(info="<html>
        <a href=\"../Documentation/HydrogenLibaryDocumnetation.pdf\">PhD project by Erasmus Rothuizen</a><br><br>
        </html>"),
          Diagram(graphics={                    Text(
                extent={{-100,112},{100,72}},
                lineColor={0,0,255},
                fillColor={255,153,0},
                fillPattern=FillPattern.Solid,
                textString=
                     "%name"),     Ellipse(
                extent={{-80,80},{80,-80}},
                lineColor={0,0,0},
                lineThickness=0.5,
                fillColor={170,255,255},
                fillPattern=FillPattern.Solid)}));
      end PressurePort;

      connector HeatFlow
        "Heat flow port passing on temperature, heat flow, pressure and area"
        import SI = Modelica.SIunits;
        SI.Temperature T "Temperature";
        flow SI.HeatFlowRate Q "Heat flow";
        SI.Pressure P "Pressure";
        Real Counter "Count pieces of walls";

        annotation (Icon(graphics={Ellipse(
                extent={{-80,80},{80,-80}},
                lineColor={0,0,0},
                lineThickness=0.5,
                fillColor={255,0,0},
                fillPattern=FillPattern.Solid)}),
                               Documentation(info="<html>
        <a href=\"../Documentation/HydrogenLibaryDocumnetation.pdf\">PhD project by Erasmus Rothuizen</a><br><br>
        </html>"),
          Diagram(graphics={                    Text(
                extent={{-104,112},{96,72}},
                lineColor={0,0,255},
                fillColor={255,153,0},
                fillPattern=FillPattern.Solid,
                textString=
                     "%name"),     Ellipse(
                extent={{-80,78},{80,-82}},
                lineColor={0,0,0},
                lineThickness=0.5,
                fillColor={255,0,0},
                fillPattern=FillPattern.Solid)}));
      end HeatFlow;

      connector HeatFlow2
        "Heat flow port passing on temperature, heat flow, pressure and area"
        import SI = Modelica.SIunits;
        SI.Temperature T "Temperature";
        flow SI.HeatFlowRate Q "Heat flow";
        SI.Pressure P "Pressure";
        Real Counter "Count pieces of walls";
        SI.MassFlowRate m_flow;
        annotation (Icon(graphics={Ellipse(
                extent={{-80,80},{80,-80}},
                lineColor={0,0,0},
                lineThickness=0.5,
                fillColor={255,255,0},
                fillPattern=FillPattern.Solid)}),
                               Documentation(info="<html>
        <a href=\"../Documentation/HydrogenLibaryDocumnetation.pdf\">PhD project by Erasmus Rothuizen</a><br><br>
        </html>"),
          Diagram(graphics={                    Text(
                extent={{-104,112},{96,72}},
                lineColor={0,0,255},
                fillColor={255,153,0},
                fillPattern=FillPattern.Solid,
                textString=
                     "%name"),     Ellipse(
                extent={{-80,78},{80,-82}},
                lineColor={0,0,0},
                lineThickness=0.5,
                fillColor={255,255,0},
                fillPattern=FillPattern.Solid)}));
      end HeatFlow2;

      connector TemperaturePort
        import SI = Modelica.SIunits;
        SI.Temperature T;

        annotation (Icon(graphics={Ellipse(
                extent={{-70,90},{90,-70}},
                lineColor={0,0,0},
                lineThickness=0.5,
                fillColor={255,0,0},
                fillPattern=FillPattern.Solid)}));
      end TemperaturePort;

      connector HeatFlowTube
        import SI = Modelica.SIunits;
       SI.Temperature T;
       flow SI.MassFlowRate m_flow;
       SI.SpecificHeatCapacity cp;
       SI.CoefficientOfHeatTransfer h;
        annotation (Icon(graphics={Ellipse(
                extent={{-70,90},{90,-70}},
                lineColor={0,0,0},
                lineThickness=0.5,
                fillColor={255,0,0},
                fillPattern=FillPattern.Solid)}));
      end HeatFlowTube;
    end Ports;

    package Sources

      model H2Source
        import SI = Modelica.SIunits;

       //replaceable package Medium = CoolProp2Modelica.Media.Hydrogen (onePhase=true)
        // constrainedby Modelica.Media.Interfaces.PartialMedium annotation (choicesAllMatching=true);
       replaceable package Medium =
          ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
          Modelica.Media.Interfaces.PartialMedium;

      Medium.ThermodynamicState state1;
      // parameter Real offset=2e6;
      outer parameter SI.Pressure P_start;
      outer SI.Pressure   FP;
      outer Real APRR;

      parameter SI.Time startTime=0;
      //SI.Time duration;
      //SI.Pressure height;

       SI.Pressure p( start=P_start);
      //SI.MassFlowRate m_flow;
      parameter SI.Temperature T = 293.15; // -60 to have -40 C
      SI.SpecificEnthalpy h;
      //SI.Pressure dp( start=0);
      Integer z(start=0);
        Ports.FlowPort portA annotation (Placement(transformation(extent={{-4,-30},
                  {12,-14}}), iconTransformation(extent={{-4,-30},{12,-14}})));

      equation
        //p=50000*(time+1);
      //height = FP;
      // duration = height/(APRR/60);
        // p = P_start + (if time < startTime then 0 else if time < (startTime + duration) then (time - startTime)*height/duration else FP);

      //Deciding APRR for controls and mass balances
      if z==0 then
        der(portA.p) = APRR/60;
      else
      portA.m_flow=0;
      end if;

      when portA.p>= FP then
        z=1;
      end when;

      portA.p=p;
        state1=Medium.setState_pT(p=p,T=T);
        h=Medium.specificEnthalpy(state1);
        portA.h_outflow=h;
      //portA.m_flow=m_flow;

        annotation (Diagram(graphics), Icon(graphics={Bitmap(extent={{-86,-48},
                    {90,118}}, fileName=
                    "modelica://HySDeP/Graphics/Hydrogen refueling station - H2 source for PCM tank.png")}));
      end H2Source;
    end Sources;

    package Properties

    end Properties;

    package HeatTransfer

      package WallPieces

        model InnerWallCell

          import SI = Modelica.SIunits;
        /*********************** Thermodynamic property call ***********************************/

                 /*    replaceable package Medium = CoolProp2Modelica.Media.Hydrogen (
        onePhase=true) constrainedby Modelica.Media.Interfaces.PartialMedium 
    annotation (choicesAllMatching=true);
*/

         replaceable package Medium =
            ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
            Modelica.Media.Interfaces.PartialMedium;
            Medium.ThermodynamicState medium;

        /******************** Connectors *****************************/
          Ports.HeatFlow2 portA "Heat flow is pos when added to the wall"
            annotation (Placement(transformation(extent={{-10,50},{10,70}},
                  rotation=0)));
          Ports.HeatFlow portB "Heat flow is pos when subtracted from the wall"
            annotation (Placement(transformation(extent={{-10,-70},{10,-50}},
                  rotation=0)),HideResult=true);

        /****************** variables *******************/
          SI.Temperature T "Temperature of the wall element";
          SI.Pressure p "Pressure of the hydrogen";
          SI.ThermalResistance R "Thermal resistance of liner";
          SI.CoefficientOfHeatTransfer h "Heat transfer coefficient";
          SI.ThermalConductivity k "Thermal conductivity of liner";
          SI.SpecificHeatCapacity cp "Specific heat capacity of liner";
          SI.Density rho "Density of liner";
          Real Tau "Dimenasionless time";
          Real Ra "Dimensionless hydrogen properties number";
          Real beta "Thermal expansion coefficient";
          Real v "kinematic viscosity";
          Real a "thermal diffusivity of hydrogen";
          Real Nu "Dimensionless heat transfer number";

        Real h_TEST;
        /****************** General parameters *******************/
          constant Real g=9.82 "acceleration due to gravity";
          SI.Length dx=x1/2;

              outer SI.Length d "inside diameter of cylinder";
              outer SI.Temperature T_amb "Ambient temperature";
              outer Real y1;
              outer SI.CoefficientOfHeatTransfer  h_i;
              outer SI.Area A;
              outer SI.Length x1;

              /*
  constant Real g=9.82 "acceleration due to gravity";
  SI.Length dx=x1/2;

       SI.Length d=0.35 "inside diameter of cylinder";
       SI.Temperature T_amb=293.15 "Ambient temperature";
       Real y1;
       SI.CoefficientOfHeatTransfer  h_i=50;
       SI.Area A=2*Modelica.Constants.pi*d/2*L+2*Modelica.Constants.pi*(d/2)^2;
       parameter SI.Length L = 0.8;
       SI.Length x1=0.03;
      
      */
        /****************** Tables *******************/
          Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            tableOnFile=true,
            table=[1,850,15,481; 2,2700,236,900; 3,1286,1.17,1578; 4,1374,1.14,1075],
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-80,0},{-60,20}})));

        /****************** Initial equations *******************/
        initial equation
          T=T_amb;

        /****************** equations *******************/
        equation
        medium=Medium.setState_pT(p, T);
        tank_prop.u=y1;
        tank_prop.y[1]=rho;
        tank_prop.y[2]=k;
        tank_prop.y[3]=cp;

         Medium.thermalConductivity(medium)=a;
         Medium.dynamicViscosity(medium)/Medium.density(medium)=v;

        //calculation of rayleighs number taken from 'Natural convection cooling of
        //rectangular and cylindrical containers' by Wenxian Lin, S.W. Armfield
        Ra=g*beta*d^3*Medium.specificHeatCapacityCp(medium)*
        Medium.density(medium)^2*abs(portA.T-T)/(v*k);
        beta=1/portA.T;
        Tau=time/(d^2/(a*Ra^(1/2)));
        Nu=0.104*Ra^(0.352);

        h_TEST = Nu*k/d;
        /*
if portA.m_flow < 0.00 then
  h=Nu*k/d;
elseif portA.m_flow > 0.00 then
  h=h_i;
else
  h=50;
end if;
*/

        if portA.m_flow > 0 then
        h = 150; // equal to "h_charging" in "HeatTransferTank"
        else

          h = 50;

        end if;

        R = dx / (A*k);

          portA.Q = (portA.T-T) / (1/(h*A));
          portB.Q = (portB.T-T) / (R/2);

          A*dx*rho*cp*der(T) = portA.Q + portB.Q;

        portA.P=p;
        portA.P=portB.P;
        portA.Counter+1=portB.Counter;
          annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}),       graphics), Icon(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics={Rectangle(extent={{-90,60},{90,-60}}, lineColor={0,
                      0,255}), Text(
                  extent={{-86,16},{86,-4}},
                  lineColor={0,0,255},
                  textString="Inner_wall_discharging")}));
        end InnerWallCell;

        model OuterWallCell "Outer wall with natural convection given by h_o"
          import SI = Modelica.SIunits;

        /******************** Connectors *****************************/
         Ports.HeatFlow portA "Heat flow is pos when added to the wall"
           annotation (Placement(transformation(extent={{-10,50},{10,70}},
                 rotation=0)));

        /****************** variables *******************/
         SI.Temperature T "Temperature of the wall element";
         SI.ThermalResistance R;
         SI.ThermalConductivity k;
         SI.SpecificHeatCapacity cp;
         SI.Density rho;
         SI.Area  A "heat flow surface area, dz*dy";
         SI.Length dx=x2/2;
         SI.HeatFlowRate Q;
          SI.HeatFlowRate Q_net;

        outer parameter Boolean Adiabatic_Wall;
        /****************** General parameters *******************/
              outer parameter SI.CoefficientOfHeatTransfer h_o;
              outer SI.Temperature T_amb;
              outer Real y1;
              outer Real x2;
              outer SI.Length d;
              outer SI.Length xLiner;
              outer SI.Length  xCFRP;
              outer SI.Length L;

        /****************** Tables *******************/
         Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            tableOnFile=true,
            table=[1,850,15,481; 2,2700,236,900; 3,1286,1.17,1578; 4,1374,1.14,1075],
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
        /****************** Initial equations *******************/
        initial equation
         T=T_amb;
        /****************** equations *******************/
        equation
        A=(d/2+xLiner+xCFRP)*2*Modelica.Constants.pi*L+2*
        (d/2+xLiner+xCFRP)^2*Modelica.Constants.pi;

        tank_prop.u=y1;
        tank_prop.y[1]=rho;
        tank_prop.y[2]=k;
        tank_prop.y[3]=cp;

         R = dx / (A*k);

        portA.Q = (portA.T-T) / (R/2);
        //Q=h_o*(T_amb-T)*A;

         A*dx*rho*cp*der(T) = portA.Q + Q;

        Q_net = portA.Q + Q;

        if Adiabatic_Wall == true then
          Q = 0;

        else

        Q=h_o*(T_amb-T)*A;

        end if;

         annotation (Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},
                    {100,100}}),      graphics), Icon(coordinateSystem(
                 preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
               graphics={Rectangle(extent={{-90,60},{90,-60}}, lineColor={0,
                     0,255}), Text(
                 extent={{-78,22},{78,-14}},
                 lineColor={0,0,255},
                  textString="Outer_Wall")}));
        end OuterWallCell;

        model LinerCell
          import SI = Modelica.SIunits;

          /******************** Connectors *****************************/
         Ports.HeatFlow portA "Heat flow is pos when added to the wall"
            annotation (Placement(transformation(extent={{-10,50},{10,70}},
                  rotation=0)));
          Ports.HeatFlow portB "Heat flow is pos when subtracted from the wall"
            annotation (Placement(transformation(extent={{-10,-70},{10,-50}},
                  rotation=0)),HideResult=true);

        /****************** variables *******************/
          SI.Temperature T "Temperature of the wall element";
          SI.ThermalResistance R;
          SI.ThermalConductivity k;
          SI.SpecificHeatCapacity cp;
          SI.Density rho;
          SI.Area A "heat flow surface area, dz*dy";
          SI.Length dx=x1;
        SI.Length Radius_local;
        /****************** General parameters *******************/
         outer SI.Temperature T_amb;
         outer Real x1;
         outer Real y1;
         outer Real t1;
         outer SI.Length d;
         outer SI.Length L;

        /****************** Tables *******************/
        Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            table=[1,850,15,481; 2,2700,167,1106; 3,1286,1.17,1578; 4,1374,1.14,
                1075],
            tableOnFile=true,
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-74,46},{-54,66}})));

        /****************** Initial equations *******************/
        initial equation
          T=T_amb;
        /****************** equations *******************/
        equation
        A=((portA.Counter-0.5)*x1+d/2)*Modelica.Constants.pi*2*L+2*
        ((portA.Counter-0.5)*x1+d/2)^2*Modelica.Constants.pi;

        Radius_local = (portA.Counter-0.5)*x1+d/2;

        tank_prop.u=y1;
        tank_prop.y[1]=rho;
        tank_prop.y[2]=k;
        tank_prop.y[3]=cp;

          R = dx / (A*k);
          portA.Q = (portA.T-T) / (R/2);
          portB.Q = (portB.T-T) / (R/2);

          A*dx*rho*cp*der(T) = portA.Q + portB.Q;

        if portA.Counter>=t1+0.5 then
          portB.Counter =0;
          else
        portA.Counter+1=portB.Counter;
        end if;

        portA.P=portB.P;
          annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}),       graphics), Icon(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics={Rectangle(extent={{-90,60},{90,-60}}, lineColor={0,
                      0,255}), Text(
                  extent={{-82,36},{84,-38}},
                  lineColor={0,0,255},
                  textString="Liner")}));
        end LinerCell;

        model TankCell
          import SI = Modelica.SIunits;

          /******************** Connectors *****************************/
         Ports.HeatFlow portA "Heat flow is pos when added to the wall"
            annotation (Placement(transformation(extent={{-10,50},{10,70}},
                  rotation=0)));
         Ports.HeatFlow portB "Heat flow is pos when subtracted from the wall"
            annotation (Placement(transformation(extent={{-10,-70},{10,-50}},
                  rotation=0)),HideResult=true);

          /****************** variables *******************/
          SI.Temperature T "Temperature of the wall element";
          SI.ThermalResistance R;
          SI.ThermalConductivity k;
          SI.SpecificHeatCapacity cp;
          SI.Density rho;

          SI.Area A "heat flow surface area, dz*dy";
          SI.Area A1;
          SI.Length dx=x2;
          SI.Length Local_Radius;
          SI.Length Local_Radius1;

          /****************** General parameters *******************/
          outer SI.Temperature T_amb;
          outer SI.Length x2;
          outer SI.Length x1;
          outer Real y2;
          outer SI.Length d;
          outer SI.Length L;
          outer Real t1;
          outer Real t2;

        /****************** Tables *******************/
        Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
            table=[1,850,15,481; 2,2700,236,900; 3,1286,1.17,1578; 4,1500,0.5,940],
            tableOnFile=true,
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-76,46},{-56,66}})));
          /****************** Initial equations *******************/
        initial equation
          T=T_amb;
          /****************** equations *******************/
        equation
        A=((portA.Counter-0.5)*x2+x1*t1+d/2)*Modelica.Constants.pi*2*L+2*
        ((portA.Counter-0.5)*x2+x1*t1+d/2)^2*Modelica.Constants.pi;

        A1 = (x2+x1*t1+d/2)*Modelica.Constants.pi*2*L+2*
        (x2+x1*t1+d/2)^2*Modelica.Constants.pi;

        Local_Radius = ((portA.Counter-0.5)*x2+x1*t1+d/2);

        Local_Radius1 = (x2+x1*t1+d/2);

        tank_prop.u=y2;
        tank_prop.y[1]=rho;
        tank_prop.y[2]=k;
        tank_prop.y[3]=cp;

          R = dx / (A*k);
          portA.Q = (portA.T-T) / (R/2);
          portB.Q = (portB.T-T) / (R/2);

        A*dx*rho*cp*der(T) = portA.Q + portB.Q;
        portA.P=portB.P;
        if portA.Counter>=t2-0.5 then
          portB.Counter =0;
          else
        portA.Counter+1=portB.Counter;
        end if;
          annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}),       graphics), Icon(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics={Rectangle(extent={{-90,60},{90,-60}}, lineColor={0,
                      0,255}), Text(
                  extent={{-68,26},{66,-32}},
                  lineColor={0,0,255},
                  textString="CFRP")}));
        end TankCell;

        model TubeCell
          import SI = Modelica.SIunits;
          import C = Modelica.Constants;

        /******************** Connectors *****************************/
            Ports.TemperaturePort                     HT1
            annotation (Placement(transformation(extent={{-56,60},{-36,80}}),
                iconTransformation(extent={{-56,60},{-36,80}})));
          Ports.TemperaturePort                     HT2
            annotation (Placement(transformation(extent={{40,60},{60,80}})));

        /****************** General parameters *******************/
         outer SI.MassFlowRate m_flow;
         outer SI.SpecificHeatCapacity cp_h2;
         outer SI.CoefficientOfHeatTransfer h;
         outer SI.CoefficientOfHeatTransfer h_amb;
         outer parameter SI.Length d_i;
         outer parameter SI.Length d_o;
         outer parameter SI.Temperature T_amb;
         outer SI.Length L;
        parameter SI.Length dx=(d_o-d_i)/12;

        /****************** variables *******************/
        flow SI.HeatFlowRate[8] Q;
        SI.Temperature[7] T;
        Real R;
        SI.SpecificHeatCapacity cp;
        SI.Density rho;
        SI.Conductivity k;

        /****************** Tables *******************/
        Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            table=[1,850,15,481; 2,2700,167,1106; 3,1286,1.17,1578; 4,1374,1.14,
                1075],
            tableOnFile=true,
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-72,-10},{-52,10}})));
        /****************** Initial equations *******************/
        initial equation
          T[1]=T_amb;
          T[3]=T_amb;
          T[5]=T_amb;
          T[7]=T_amb;

        /****************** equations *******************/
        equation

          tank_prop.u=1;
          tank_prop.y[1]=rho;
          tank_prop.y[2]=k;
          tank_prop.y[3]=cp;

          R=dx/((d_i*C.pi*L)*k);
        //Between hydrogen and inner wall
          Q[1]=(HT1.T-T[1])/(1/(h*(d_i*C.pi*L)));
          Q[2]=(T[2]-T[1])/(R/2);
          (d_i*C.pi*L)*dx*rho*cp*der(T[1]) = Q[1] + Q[2];
        //Between two wall volumes
          -Q[2]=Q[3];
          Q[3]=(T[2]-T[3])/(R/2);
          Q[4]=(T[4]-T[3])/(R/2);
          ((d_i/2+dx*2)*2*C.pi*L)*2*dx*rho*cp*der(T[3]) = Q[3] + Q[4];
        //Between two wall volumes
          -Q[4]=Q[5];
          Q[5]=(T[4]-T[5])/(R/2);
          Q[6]=(T[6]-T[5])/(R/2);
          ((d_i/2+dx*4)*2*C.pi*L)*2*dx*rho*cp*der(T[5]) = Q[5] + Q[6];
        //Between outer wall volume and ambient
          -Q[6]=Q[7];
          Q[7]=(T[6]-T[7])/(R/2);
          Q[8]=h_amb*(T_amb-T[7])*((d_i/2+dx*6)*2*C.pi*L);
          ((d_i/2+dx*6)*2*C.pi*L)*dx*rho*cp*der(T[7]) = Q[7] + Q[8];

          -Q[1]=cp_h2*m_flow*(HT1.T-HT2.T);

          annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                    -100},{100,100}}), graphics), Icon(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics={Rectangle(extent={{-84,70},{96,-50}}, lineColor={0,
                      0,255}), Text(
                  extent={{-56,24},{66,-28}},
                  lineColor={0,0,255},
                  textString="Tube cell
")}));
        end TubeCell;

        model Liner5Pieces
        //Import of models
          LinerCell
                liner
            annotation (Placement(transformation(extent={{-10,50},{10,70}})));
          LinerCell
                liner1
            annotation (Placement(transformation(extent={{-10,24},{10,44}})));
          LinerCell
                liner2
            annotation (Placement(transformation(extent={{-10,0},{10,20}})));
          LinerCell
                liner3
            annotation (Placement(transformation(extent={{-10,-24},{10,-4}})));
          LinerCell
                liner4
            annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
          Ports.HeatFlow heatFlow
            annotation (Placement(transformation(extent={{-10,76},{10,96}})));
         Ports.HeatFlow heatFlow1
            annotation (Placement(transformation(extent={{-10,-80},{10,-60}})));
        equation
          connect(liner.portB, liner1.portA) annotation (Line(
              points={{0,54},{0,40}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(liner1.portB, liner2.portA) annotation (Line(
              points={{0,28},{0,16}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(liner2.portB, liner3.portA) annotation (Line(
              points={{0,4},{0,-8}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(liner3.portB, liner4.portA) annotation (Line(
              points={{0,-20},{0,-34}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(liner.portA, heatFlow) annotation (Line(
              points={{0,66},{0,86}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(liner4.portB, heatFlow1) annotation (Line(
              points={{0,-46},{0,-70}},
              color={0,0,255},
              smooth=Smooth.None));
          annotation (Diagram(graphics));
        end Liner5Pieces;

        model Tank10Pieces
        //import models
          TankCell
               cFRP
            annotation (Placement(transformation(extent={{-50,36},{-30,56}})));
          TankCell
               cFRP1
            annotation (Placement(transformation(extent={{-50,12},{-30,32}})));
          TankCell
               cFRP2
            annotation (Placement(transformation(extent={{-50,-12},{-30,8}})));
          TankCell
               cFRP3
            annotation (Placement(transformation(extent={{-50,-36},{-30,-16}})));
          TankCell
               cFRP4
            annotation (Placement(transformation(extent={{-50,-60},{-30,-40}})));
          TankCell
               cFRP5
            annotation (Placement(transformation(extent={{-10,36},{10,56}})));
          TankCell
               cFRP6
            annotation (Placement(transformation(extent={{-10,12},{10,32}})));
          TankCell
               cFRP7
            annotation (Placement(transformation(extent={{-10,-12},{10,8}})));
          TankCell
               cFRP8
            annotation (Placement(transformation(extent={{-10,-36},{10,-16}})));
          TankCell
               cFRP9
            annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));
          Ports.HeatFlow heatFlow
            annotation (Placement(transformation(extent={{-10,74},{10,94}})));
          Ports.HeatFlow heatFlow1
            annotation (Placement(transformation(extent={{-10,-88},{10,-68}})));
        equation
          connect(cFRP.portA, heatFlow) annotation (Line(
              points={{-40,52},{-40,84},{0,84}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP.portB, cFRP1.portA) annotation (Line(
              points={{-40,40},{-40,28}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP1.portB, cFRP2.portA) annotation (Line(
              points={{-40,16},{-40,4}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP2.portB, cFRP3.portA) annotation (Line(
              points={{-40,-8},{-40,-20}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP3.portB, cFRP4.portA) annotation (Line(
              points={{-40,-32},{-40,-44}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP4.portB, cFRP5.portA) annotation (Line(
              points={{-40,-56},{-40,-64},{-18,-64},{-18,60},{0,60},{0,52}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP5.portB, cFRP6.portA) annotation (Line(
              points={{0,40},{0,28}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP6.portB, cFRP7.portA) annotation (Line(
              points={{0,16},{0,4}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP7.portB, cFRP8.portA) annotation (Line(
              points={{0,-8},{0,-20}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP8.portB, cFRP9.portA) annotation (Line(
              points={{0,-32},{0,-44}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP9.portB, heatFlow1) annotation (Line(
              points={{0,-56},{0,-56},{0,-78}},
              color={0,0,255},
              smooth=Smooth.None));
          annotation (Diagram(graphics));
        end Tank10Pieces;

        model Tube5Pieces

          TubeCell   tubeHeatTransfer3Testing
            annotation (Placement(transformation(extent={{-60,4},{-40,24}})));
          TubeCell   tubeHeatTransfer3Testing1
            annotation (Placement(transformation(extent={{-34,4},{-14,24}})));
          TubeCell   tubeHeatTransfer3Testing2
            annotation (Placement(transformation(extent={{-8,4},{12,24}})));
          TubeCell  tubeHeatTransfer3Testing3
            annotation (Placement(transformation(extent={{18,4},{38,24}})));
          TubeCell  tubeHeatTransfer3Testing4
            annotation (Placement(transformation(extent={{44,4},{64,24}})));

          Ports.TemperaturePort                     HT1
            annotation (Placement(transformation(extent={{-72,50},{-52,70}}),
                iconTransformation(extent={{-72,50},{-52,70}})));

          Ports.TemperaturePort                     HT2
            annotation (Placement(transformation(extent={{48,50},{68,70}}),
                iconTransformation(extent={{48,50},{68,70}})));
        equation
          connect(HT1,tubeHeatTransfer3Testing. HT1) annotation (Line(
              points={{-62,60},{-54.6,60},{-54.6,21}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tubeHeatTransfer3Testing1.HT2,tubeHeatTransfer3Testing2. HT1)
            annotation (Line(
              points={{-19,21},{-10,22},{-2.6,21}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tubeHeatTransfer3Testing2.HT2,tubeHeatTransfer3Testing3. HT1)
            annotation (Line(
              points={{7,21},{7,18.5},{23.4,18.5},{23.4,21}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tubeHeatTransfer3Testing3.HT2,tubeHeatTransfer3Testing4. HT1)
            annotation (Line(
              points={{33,21},{49.4,21}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tubeHeatTransfer3Testing4.HT2, HT2) annotation (Line(
              points={{59,21},{59,28.5},{58,28.5},{58,60}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tubeHeatTransfer3Testing1.HT1, tubeHeatTransfer3Testing.HT2)
            annotation (Line(
              points={{-28.6,21},{-45,21}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                    -100},{100,100}}), graphics), Icon(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics));
        end Tube5Pieces;
      end WallPieces;

      package WallPieces_Area

        model InnerWallCell

          import SI = Modelica.SIunits;
        /*********************** Thermodynamic property call ***********************************/

                 /*    replaceable package Medium = CoolProp2Modelica.Media.Hydrogen (
        onePhase=true) constrainedby Modelica.Media.Interfaces.PartialMedium 
    annotation (choicesAllMatching=true);
*/

         replaceable package Medium =
            ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
            Modelica.Media.Interfaces.PartialMedium;
            Medium.ThermodynamicState medium;

        /******************** Connectors *****************************/
          Ports.HeatFlow2 portA "Heat flow is pos when added to the wall"
            annotation (Placement(transformation(extent={{-10,50},{10,70}},
                  rotation=0)));
          Ports.HeatFlow portB "Heat flow is pos when subtracted from the wall"
            annotation (Placement(transformation(extent={{-10,-70},{10,-50}},
                  rotation=0)),HideResult=true);

        /****************** variables *******************/
          SI.Temperature T "Temperature of the wall element";
          SI.Pressure p "Pressure of the hydrogen";
          SI.ThermalResistance R "Thermal resistance of liner";
          SI.CoefficientOfHeatTransfer h "Heat transfer coefficient";
          SI.ThermalConductivity k "Thermal conductivity of liner";
          SI.SpecificHeatCapacity cp "Specific heat capacity of liner";
          SI.Density rho "Density of liner";
          Real Tau "Dimenasionless time";
          Real Ra "Dimensionless hydrogen properties number";
          Real beta "Thermal expansion coefficient";
          Real v "kinematic viscosity";
          Real a "thermal diffusivity of hydrogen";
          Real Nu "Dimensionless heat transfer number";

        /****************** General parameters *******************/
          constant Real g=9.82 "acceleration due to gravity";
          SI.Length dx=x1/2;

              outer SI.Length d "inside diameter of cylinder";
              outer SI.Temperature T_amb "Ambient temperature";
              outer Real y1;
              outer SI.CoefficientOfHeatTransfer  h_i;
              outer SI.Area A;
              outer SI.Length x1;

              /*
  constant Real g=9.82 "acceleration due to gravity";
  SI.Length dx=x1/2;

       SI.Length d=0.35 "inside diameter of cylinder";
       SI.Temperature T_amb=293.15 "Ambient temperature";
       Real y1;
       SI.CoefficientOfHeatTransfer  h_i=50;
       SI.Area A=2*Modelica.Constants.pi*d/2*L+2*Modelica.Constants.pi*(d/2)^2;
       parameter SI.Length L = 0.8;
       SI.Length x1=0.03;
      
      */
        /****************** Tables *******************/
          Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            tableOnFile=true,
            table=[1,850,15,481; 2,2700,236,900; 3,1286,1.17,1578; 4,1374,1.14,1075],
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-80,0},{-60,20}})));

        /****************** Initial equations *******************/
        initial equation
          T=T_amb;

        /****************** equations *******************/
        equation
        medium=Medium.setState_pT(p, T);
        tank_prop.u=y1;
        tank_prop.y[1]=rho;
        tank_prop.y[2]=k;
        tank_prop.y[3]=cp;

         Medium.thermalConductivity(medium)=a;
         Medium.dynamicViscosity(medium)/Medium.density(medium)=v;

        //calculation of rayleighs number taken from 'Natural convection cooling of
        //rectangular and cylindrical containers' by Wenxian Lin, S.W. Armfield
        Ra=g*beta*d^3*Medium.specificHeatCapacityCp(medium)*
        Medium.density(medium)^2*abs(portA.T-T)/(v*k);
        beta=1/portA.T;
        Tau=time/(d^2/(a*Ra^(1/2)));
        Nu=0.104*Ra^(0.352);

        if portA.m_flow < 0.00 then
          h=Nu*k/d;
        elseif portA.m_flow > 0.00 then
          h=h_i;
        else
          h=50;
        end if;

        R = dx / (A*k);

          portA.Q = (portA.T-T) / (1/(h*A));
          portB.Q = (portB.T-T) / R; // instead of R/2, because there is only one conduction resistance and dx = x1/2 already

          A*dx*rho*cp*der(T) = portA.Q + portB.Q;

        portA.P=p;
        portA.P=portB.P;
        portA.Counter+1=portB.Counter;
          annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}),       graphics), Icon(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics={Rectangle(extent={{-90,60},{90,-60}}, lineColor={0,
                      0,255}), Text(
                  extent={{-86,16},{86,-4}},
                  lineColor={0,0,255},
                  textString="Inner_wall_discharging")}));
        end InnerWallCell;

        model OuterWallCell "Outer wall with natural convection given by h_o"
          import SI = Modelica.SIunits;

        /******************** Connectors *****************************/
         Ports.HeatFlow portA "Heat flow is pos when added to the wall"
           annotation (Placement(transformation(extent={{-10,50},{10,70}},
                 rotation=0)));

        /****************** variables *******************/
         SI.Temperature T "Temperature of the wall element";
         SI.ThermalResistance R;
         SI.ThermalConductivity k;
         SI.SpecificHeatCapacity cp;
         SI.Density rho;
         SI.Area  A "heat flow surface area, dz*dy";
         SI.Length dx=x2/2;
         SI.HeatFlowRate Q;
          SI.HeatFlowRate Q_net;
        SI.Length Local_Radius;
        outer parameter Boolean Adiabatic_Wall;
        /****************** General parameters *******************/
              outer parameter SI.CoefficientOfHeatTransfer h_o;
              outer SI.Temperature T_amb;
              outer Real y1;
              outer Real x2;
              outer SI.Length d;
              outer SI.Length xLiner;
              outer SI.Length  xCFRP;
              outer SI.Length L;

        /****************** Tables *******************/
         Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            tableOnFile=true,
            table=[1,850,15,481; 2,2700,236,900; 3,1286,1.17,1578; 4,1374,1.14,1075],
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
        /****************** Initial equations *******************/
        initial equation
         T=T_amb;
        /****************** equations *******************/
        equation
        A=(d/2+xLiner+xCFRP+xLiner/10)*2*Modelica.Constants.pi*L+2*
        (d/2+xLiner+xCFRP)^2*Modelica.Constants.pi;

        Local_Radius = (d/2+xLiner+xCFRP+xLiner/10); // xLiner/10 accounts for the half thermal resistance in "InnerWallcell"

        tank_prop.u=y1;
        tank_prop.y[1]=rho;
        tank_prop.y[2]=k;
        tank_prop.y[3]=cp;

         R = dx / (A*k);

        portA.Q = (portA.T-T) / (R/2);
        //Q=h_o*(T_amb-T)*A;

         A*dx*rho*cp*der(T) = portA.Q + Q;

        Q_net = portA.Q + Q;

        if Adiabatic_Wall == true then
          Q = 0;

        else

        Q=h_o*(T_amb-T)*A;

        end if;

         annotation (Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},
                    {100,100}}),      graphics), Icon(coordinateSystem(
                 preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
               graphics={Rectangle(extent={{-90,60},{90,-60}}, lineColor={0,
                     0,255}), Text(
                 extent={{-78,22},{78,-14}},
                 lineColor={0,0,255},
                  textString="Outer_Wall")}));
        end OuterWallCell;

        model LinerCell
          import SI = Modelica.SIunits;

          /******************** Connectors *****************************/
         Ports.HeatFlow portA "Heat flow is pos when added to the wall"
            annotation (Placement(transformation(extent={{-10,50},{10,70}},
                  rotation=0)));
          Ports.HeatFlow portB "Heat flow is pos when subtracted from the wall"
            annotation (Placement(transformation(extent={{-10,-70},{10,-50}},
                  rotation=0)),HideResult=true);

        /****************** variables *******************/
          SI.Temperature T "Temperature of the wall element";
          SI.ThermalResistance R;
          SI.ThermalConductivity k;
          SI.SpecificHeatCapacity cp;
          SI.Density rho;
          SI.Area A "heat flow surface area, dz*dy";
          SI.Length dx=x1;
        SI.Length Radius_local;
        /****************** General parameters *******************/
         outer SI.Temperature T_amb;
         outer Real x1;
         outer Real y1;
         outer Real t1;
         outer SI.Length d;
         outer SI.Length L;

        /****************** Tables *******************/
        Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            table=[1,850,15,481; 2,2700,167,1106; 3,1286,1.17,1578; 4,1374,1.14,
                1075],
            tableOnFile=true,
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-74,46},{-54,66}})));

        /****************** Initial equations *******************/
        initial equation
          T=T_amb;
        /****************** equations *******************/
        equation
        A=((portA.Counter-0.5)*x1+d/2)*Modelica.Constants.pi*2*L+2*
        ((portA.Counter-0.5)*x1+d/2)^2*Modelica.Constants.pi;

        Radius_local = (portA.Counter-0.5)*x1+d/2;

        tank_prop.u=y1;
        tank_prop.y[1]=rho;
        tank_prop.y[2]=k;
        tank_prop.y[3]=cp;

          R = dx / (A*k);
          portA.Q = (portA.T-T) / (R/2);
          portB.Q = (portB.T-T) / (R/2);

          A*dx*rho*cp*der(T) = portA.Q + portB.Q;

        if portA.Counter>=t1+0.5 then
          portB.Counter =0;
          else
        portA.Counter+1=portB.Counter;
        end if;

        portA.P=portB.P;
          annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}),       graphics), Icon(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics={Rectangle(extent={{-90,60},{90,-60}}, lineColor={0,
                      0,255}), Text(
                  extent={{-82,36},{84,-38}},
                  lineColor={0,0,255},
                  textString="Liner")}));
        end LinerCell;

        model TankCell
          import SI = Modelica.SIunits;

          /******************** Connectors *****************************/
         Ports.HeatFlow portA "Heat flow is pos when added to the wall"
            annotation (Placement(transformation(extent={{-10,50},{10,70}},
                  rotation=0)));
         Ports.HeatFlow portB "Heat flow is pos when subtracted from the wall"
            annotation (Placement(transformation(extent={{-10,-70},{10,-50}},
                  rotation=0)),HideResult=true);

          /****************** variables *******************/
          SI.Temperature T "Temperature of the wall element";
          SI.ThermalResistance R;
          SI.ThermalConductivity k;
          SI.SpecificHeatCapacity cp;
          SI.Density rho;

          SI.Area A "heat flow surface area, dz*dy";
          SI.Length dx=x2;
          SI.Length Local_Radius;

          /****************** General parameters *******************/
          outer SI.Temperature T_amb;
          outer SI.Length x2;
          outer SI.Length x1;
          outer Real y2;
          outer SI.Length d;
          outer SI.Length L;
          outer Real t1;
          outer Real t2;

        /****************** Tables *******************/
        Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
            table=[1,850,15,481; 2,2700,236,900; 3,1286,1.17,1578; 4,1500,0.5,940],
            tableOnFile=true,
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-76,46},{-56,66}})));
          /****************** Initial equations *******************/
        initial equation
          T=T_amb;
          /****************** equations *******************/
        equation
        A=((portA.Counter)*x2+x1*t1+d/2)*Modelica.Constants.pi*2*L+2*
        ((portA.Counter)*x2+x1*t1+d/2)^2*Modelica.Constants.pi;

        Local_Radius = (portA.Counter*x2) +x1*t1+d/2;

        tank_prop.u=y2;
        tank_prop.y[1]=rho;
        tank_prop.y[2]=k;
        tank_prop.y[3]=cp;

          R = dx / (A*k);
          portA.Q = (portA.T-T) / (R/2);
          portB.Q = (portB.T-T) / (R/2);

        A*dx*rho*cp*der(T) = portA.Q + portB.Q;
        portA.P=portB.P;
        if portA.Counter>=t2-0.5 then
          portB.Counter =0;
          else
        portA.Counter+1=portB.Counter;
        end if;
          annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}),       graphics), Icon(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics={Rectangle(extent={{-90,60},{90,-60}}, lineColor={0,
                      0,255}), Text(
                  extent={{-68,26},{66,-32}},
                  lineColor={0,0,255},
                  textString="CFRP")}));
        end TankCell;

        model TankCell_Inner
          import SI = Modelica.SIunits;

          /******************** Connectors *****************************/
         Ports.HeatFlow portA "Heat flow is pos when added to the wall"
            annotation (Placement(transformation(extent={{-10,50},{10,70}},
                  rotation=0)));
         Ports.HeatFlow portB "Heat flow is pos when subtracted from the wall"
            annotation (Placement(transformation(extent={{-10,-70},{10,-50}},
                  rotation=0)),HideResult=true);

          /****************** variables *******************/
          SI.Temperature T "Temperature of the wall element";
          SI.ThermalResistance R;
          SI.ThermalConductivity k;
          SI.SpecificHeatCapacity cp;
          SI.Density rho;

          SI.Area A "heat flow surface area, dz*dy";
          SI.Length dx=x2;
          SI.Length Local_Radius;

          /****************** General parameters *******************/
          outer SI.Temperature T_amb;
          outer SI.Length x2;
          outer SI.Length x1;
          outer Real y2;
          outer SI.Length d;
          outer SI.Length L;
          outer Real t1;
          outer Real t2;

        /****************** Tables *******************/
        Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
            table=[1,850,15,481; 2,2700,236,900; 3,1286,1.17,1578; 4,1500,0.5,940],
            tableOnFile=true,
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-76,46},{-56,66}})));
          /****************** Initial equations *******************/
        initial equation
          T=T_amb;

          /****************** equations *******************/
        equation
        A=(x1*t1+d/2)*Modelica.Constants.pi*2*L+2*
        (x1*t1+d/2)^2*Modelica.Constants.pi;

        Local_Radius = (x1*t1+d/2);

        tank_prop.u=y2;
        tank_prop.y[1]=rho;
        tank_prop.y[2]=k;
        tank_prop.y[3]=cp;

          R = dx / (A*k);
          portA.Q = (portA.T-T) / (R/2);
          portB.Q = (portB.T-T) / (R/2);

        A*dx*rho*cp*der(T) = portA.Q + portB.Q;
        portA.P=portB.P;

        portA.Counter-portA.Counter +1=portB.Counter; // Restart portB.Counter = 1 for CFRP

        // if portA.Counter>=t2-0.5 then
        //   portB.Counter =0;
        //   else
        // portA.Counter+1=portB.Counter;
        // end if;

          annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}),       graphics), Icon(coordinateSystem(
                  preserveAspectRatio=true,  extent={{-100,-100},{100,100}}),
                graphics={Rectangle(extent={{-90,60},{90,-60}}, lineColor={0,
                      0,255}), Text(
                  extent={{-68,26},{66,-32}},
                  lineColor={0,0,255},
                  textString="CFRP"),
                Text(
                  extent={{-88,54},{-16,20}},
                  lineColor={0,0,255},
                  textString="Inner")}));
        end TankCell_Inner;

        model TubeCell
          import SI = Modelica.SIunits;
          import C = Modelica.Constants;

        /******************** Connectors *****************************/
            Ports.TemperaturePort                     HT1
            annotation (Placement(transformation(extent={{-56,60},{-36,80}}),
                iconTransformation(extent={{-56,60},{-36,80}})));
          Ports.TemperaturePort                     HT2
            annotation (Placement(transformation(extent={{40,60},{60,80}})));

        /****************** General parameters *******************/
         outer SI.MassFlowRate m_flow;
         outer SI.SpecificHeatCapacity cp_h2;
         outer SI.CoefficientOfHeatTransfer h;
         outer SI.CoefficientOfHeatTransfer h_amb;
         outer parameter SI.Length d_i;
         outer parameter SI.Length d_o;
         outer parameter SI.Temperature T_amb;
         outer SI.Length L;
        parameter SI.Length dx=(d_o-d_i)/12;

        /****************** variables *******************/
        flow SI.HeatFlowRate[8] Q;
        SI.Temperature[7] T;
        Real R;
        SI.SpecificHeatCapacity cp;
        SI.Density rho;
        SI.Conductivity k;

        /****************** Tables *******************/
        Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            table=[1,850,15,481; 2,2700,167,1106; 3,1286,1.17,1578; 4,1374,1.14,
                1075],
            tableOnFile=true,
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-72,-10},{-52,10}})));
        /****************** Initial equations *******************/
        initial equation
          T[1]=T_amb;
          T[3]=T_amb;
          T[5]=T_amb;
          T[7]=T_amb;

        /****************** equations *******************/
        equation

          tank_prop.u=1;
          tank_prop.y[1]=rho;
          tank_prop.y[2]=k;
          tank_prop.y[3]=cp;

          R=dx/((d_i*C.pi*L)*k);
        //Between hydrogen and inner wall
          Q[1]=(HT1.T-T[1])/(1/(h*(d_i*C.pi*L)));
          Q[2]=(T[2]-T[1])/(R/2);
          (d_i*C.pi*L)*dx*rho*cp*der(T[1]) = Q[1] + Q[2];
        //Between two wall volumes
          -Q[2]=Q[3];
          Q[3]=(T[2]-T[3])/(R/2);
          Q[4]=(T[4]-T[3])/(R/2);
          ((d_i/2+dx*2)*2*C.pi*L)*2*dx*rho*cp*der(T[3]) = Q[3] + Q[4];
        //Between two wall volumes
          -Q[4]=Q[5];
          Q[5]=(T[4]-T[5])/(R/2);
          Q[6]=(T[6]-T[5])/(R/2);
          ((d_i/2+dx*4)*2*C.pi*L)*2*dx*rho*cp*der(T[5]) = Q[5] + Q[6];
        //Between outer wall volume and ambient
          -Q[6]=Q[7];
          Q[7]=(T[6]-T[7])/(R/2);
          Q[8]=h_amb*(T_amb-T[7])*((d_i/2+dx*6)*2*C.pi*L);
          ((d_i/2+dx*6)*2*C.pi*L)*dx*rho*cp*der(T[7]) = Q[7] + Q[8];

          -Q[1]=cp_h2*m_flow*(HT1.T-HT2.T);

          annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                    -100},{100,100}}), graphics), Icon(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics={Rectangle(extent={{-84,70},{96,-50}}, lineColor={0,
                      0,255}), Text(
                  extent={{-56,24},{66,-28}},
                  lineColor={0,0,255},
                  textString="Tube cell
")}));
        end TubeCell;

        model Liner5Pieces
        //Import of models
          LinerCell
                liner
            annotation (Placement(transformation(extent={{-10,50},{10,70}})));
          LinerCell
                liner1
            annotation (Placement(transformation(extent={{-10,24},{10,44}})));
          LinerCell
                liner2
            annotation (Placement(transformation(extent={{-10,0},{10,20}})));
          LinerCell
                liner3
            annotation (Placement(transformation(extent={{-10,-24},{10,-4}})));
          LinerCell
                liner4
            annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
          Ports.HeatFlow heatFlow
            annotation (Placement(transformation(extent={{-10,76},{10,96}})));
         Ports.HeatFlow heatFlow1
            annotation (Placement(transformation(extent={{-10,-80},{10,-60}})));
        equation
          connect(liner.portB, liner1.portA) annotation (Line(
              points={{0,54},{0,40}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(liner1.portB, liner2.portA) annotation (Line(
              points={{0,28},{0,16}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(liner2.portB, liner3.portA) annotation (Line(
              points={{0,4},{0,-8}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(liner3.portB, liner4.portA) annotation (Line(
              points={{0,-20},{0,-34}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(liner.portA, heatFlow) annotation (Line(
              points={{0,66},{0,86}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(liner4.portB, heatFlow1) annotation (Line(
              points={{0,-46},{0,-70}},
              color={0,0,255},
              smooth=Smooth.None));
          annotation (Diagram(graphics));
        end Liner5Pieces;

        model Tank10Pieces
        //import models
          TankCell
               cFRP1
            annotation (Placement(transformation(extent={{-50,12},{-30,32}})));
          TankCell
               cFRP2
            annotation (Placement(transformation(extent={{-50,-12},{-30,8}})));
          TankCell
               cFRP3
            annotation (Placement(transformation(extent={{-50,-36},{-30,-16}})));
          TankCell
               cFRP4
            annotation (Placement(transformation(extent={{-50,-60},{-30,-40}})));
          TankCell
               cFRP5
            annotation (Placement(transformation(extent={{-10,36},{10,56}})));
          TankCell
               cFRP6
            annotation (Placement(transformation(extent={{-10,12},{10,32}})));
          TankCell
               cFRP7
            annotation (Placement(transformation(extent={{-10,-12},{10,8}})));
          TankCell
               cFRP8
            annotation (Placement(transformation(extent={{-10,-36},{10,-16}})));
          TankCell
               cFRP9
            annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));
          Ports.HeatFlow heatFlow
            annotation (Placement(transformation(extent={{-10,74},{10,94}})));
          Ports.HeatFlow heatFlow1
            annotation (Placement(transformation(extent={{-10,-88},{10,-68}})));
          TankCell_Inner tankCell_Inner
            annotation (Placement(transformation(extent={{-50,36},{-30,56}})));
        equation
          connect(cFRP1.portB, cFRP2.portA) annotation (Line(
              points={{-40,16},{-40,4}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP2.portB, cFRP3.portA) annotation (Line(
              points={{-40,-8},{-40,-20}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP3.portB, cFRP4.portA) annotation (Line(
              points={{-40,-32},{-40,-44}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP4.portB, cFRP5.portA) annotation (Line(
              points={{-40,-56},{-40,-64},{-18,-64},{-18,60},{0,60},{0,52}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP5.portB, cFRP6.portA) annotation (Line(
              points={{0,40},{0,28}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP6.portB, cFRP7.portA) annotation (Line(
              points={{0,16},{0,4}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP7.portB, cFRP8.portA) annotation (Line(
              points={{0,-8},{0,-20}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP8.portB, cFRP9.portA) annotation (Line(
              points={{0,-32},{0,-44}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP9.portB, heatFlow1) annotation (Line(
              points={{0,-56},{0,-56},{0,-78}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(tankCell_Inner.portB, cFRP1.portA) annotation (Line(
              points={{-40,40},{-40,28}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tankCell_Inner.portA, heatFlow) annotation (Line(
              points={{-40,52},{-40,84},{0,84}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          annotation (Diagram(graphics));
        end Tank10Pieces;

        model Tube5Pieces

          TubeCell   tubeHeatTransfer3Testing
            annotation (Placement(transformation(extent={{-60,4},{-40,24}})));
          TubeCell   tubeHeatTransfer3Testing1
            annotation (Placement(transformation(extent={{-34,4},{-14,24}})));
          TubeCell   tubeHeatTransfer3Testing2
            annotation (Placement(transformation(extent={{-8,4},{12,24}})));
          TubeCell  tubeHeatTransfer3Testing3
            annotation (Placement(transformation(extent={{18,4},{38,24}})));
          TubeCell  tubeHeatTransfer3Testing4
            annotation (Placement(transformation(extent={{44,4},{64,24}})));

          Ports.TemperaturePort                     HT1
            annotation (Placement(transformation(extent={{-72,50},{-52,70}}),
                iconTransformation(extent={{-72,50},{-52,70}})));

          Ports.TemperaturePort                     HT2
            annotation (Placement(transformation(extent={{48,50},{68,70}}),
                iconTransformation(extent={{48,50},{68,70}})));
        equation
          connect(HT1,tubeHeatTransfer3Testing. HT1) annotation (Line(
              points={{-62,60},{-54.6,60},{-54.6,21}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tubeHeatTransfer3Testing1.HT2,tubeHeatTransfer3Testing2. HT1)
            annotation (Line(
              points={{-19,21},{-10,22},{-2.6,21}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tubeHeatTransfer3Testing2.HT2,tubeHeatTransfer3Testing3. HT1)
            annotation (Line(
              points={{7,21},{7,18.5},{23.4,18.5},{23.4,21}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tubeHeatTransfer3Testing3.HT2,tubeHeatTransfer3Testing4. HT1)
            annotation (Line(
              points={{33,21},{49.4,21}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tubeHeatTransfer3Testing4.HT2, HT2) annotation (Line(
              points={{59,21},{59,28.5},{58,28.5},{58,60}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tubeHeatTransfer3Testing1.HT1, tubeHeatTransfer3Testing.HT2)
            annotation (Line(
              points={{-28.6,21},{-45,21}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                    -100},{100,100}}), graphics), Icon(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics));
        end Tube5Pieces;
      end WallPieces_Area;

      package WallPieces_Area_middle

        model InnerWallCell

          import SI = Modelica.SIunits;
        /*********************** Thermodynamic property call ***********************************/

                 /*    replaceable package Medium = CoolProp2Modelica.Media.Hydrogen (
        onePhase=true) constrainedby Modelica.Media.Interfaces.PartialMedium 
    annotation (choicesAllMatching=true);
*/

         replaceable package Medium =
            ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
            Modelica.Media.Interfaces.PartialMedium;
            Medium.ThermodynamicState medium;

        /******************** Connectors *****************************/
          Ports.HeatFlow2 portA "Heat flow is pos when added to the wall"
            annotation (Placement(transformation(extent={{-10,50},{10,70}},
                  rotation=0)));
          Ports.HeatFlow portB "Heat flow is pos when subtracted from the wall"
            annotation (Placement(transformation(extent={{-10,-70},{10,-50}},
                  rotation=0)),HideResult=true);

        /****************** variables *******************/
          SI.Temperature T "Temperature of the wall element";
          SI.Pressure p "Pressure of the hydrogen";
          SI.ThermalResistance R "Thermal resistance of liner";
          SI.CoefficientOfHeatTransfer h "Heat transfer coefficient";
          SI.ThermalConductivity k "Thermal conductivity of liner";
          SI.SpecificHeatCapacity cp "Specific heat capacity of liner";
          SI.Density rho "Density of liner";
          Real Tau "Dimenasionless time";
          Real Ra "Dimensionless hydrogen properties number";
          Real beta "Thermal expansion coefficient";
          Real v "kinematic viscosity";
          Real a "thermal diffusivity of hydrogen";
          Real Nu "Dimensionless heat transfer number";

        Real h_TEST;

        /****************** General parameters *******************/
          constant Real g=9.82 "acceleration due to gravity";
          SI.Length dx=x1/2;

              outer SI.Length d "inside diameter of cylinder";
              outer SI.Temperature T_amb "Ambient temperature";
              outer Real y1;
              outer parameter SI.CoefficientOfHeatTransfer  h_i;
              outer SI.Area A;
              outer SI.Length x1;

        /****************** Tables *******************/
          Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            tableOnFile=true,
            table=[1,850,15,481; 2,2700,236,900; 3,1286,1.17,1578; 4,1374,1.14,1075],
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-80,0},{-60,20}})));

        /****************** Initial equations *******************/
        initial equation
          T=T_amb;

        /****************** equations *******************/
        equation
        medium=Medium.setState_pT(p, T);
        tank_prop.u=y1;
        tank_prop.y[1]=rho;
        tank_prop.y[2]=k;
        tank_prop.y[3]=cp;

         Medium.thermalConductivity(medium)=a;
         Medium.dynamicViscosity(medium)/Medium.density(medium)=v;

        //calculation of rayleighs number taken from 'Natural convection cooling of
        //rectangular and cylindrical containers' by Wenxian Lin, S.W. Armfield
        Ra=g*beta*d^3*Medium.specificHeatCapacityCp(medium)*
        Medium.density(medium)^2*abs(portA.T-T)/(v*k);
        beta=1/portA.T;
        Tau=time/(d^2/(a*Ra^(1/2)));
        Nu=0.104*Ra^(0.352);

        h_TEST = Nu*k/d;

        /*
if portA.m_flow < 0.00 then
  h=Nu*k/d;
elseif portA.m_flow > 0.00 then
  h=h_i;
else
  h=50;
end if;
*/

        if portA.m_flow > 0 then
        h = h_i; // equal to "h_charging" in "HeatTransferTank"
        else

          h = 50;

        end if;

        R = dx / (A*k);

          portA.Q = (portA.T-T) / (1/(h*A));
          portB.Q = (portB.T-T) / R; // instead of R/2, because there is only one conduction resistance and dx = x1/2 already

          A*dx*rho*cp*der(T) = portA.Q + portB.Q;

        portA.P=p;
        portA.P=portB.P;
        portA.Counter+1=portB.Counter;
          annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}),       graphics), Icon(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics={Rectangle(extent={{-90,60},{90,-60}}, lineColor={0,
                      0,255}), Text(
                  extent={{-86,16},{86,-4}},
                  lineColor={0,0,255},
                  textString="Inner_wall_discharging")}));
        end InnerWallCell;

        model OuterWallCell "Outer wall with natural convection given by h_o"
          import SI = Modelica.SIunits;

        /******************** Connectors *****************************/
         Ports.HeatFlow portA "Heat flow is pos when added to the wall"
           annotation (Placement(transformation(extent={{-10,50},{10,70}},
                 rotation=0)));

        /****************** variables *******************/
         SI.Temperature T "Temperature of the wall element";
         SI.ThermalResistance R;
         SI.ThermalConductivity k;
         SI.SpecificHeatCapacity cp;
         SI.Density rho;
         SI.Area  A "heat flow surface area, dz*dy";
         SI.Length dx=x2/2;
         SI.HeatFlowRate Q;
          SI.HeatFlowRate Q_net;
        SI.Length Local_Radius;
        outer parameter Boolean Adiabatic_Wall;
        /****************** General parameters *******************/
              outer parameter SI.CoefficientOfHeatTransfer h_o;
              outer SI.Temperature T_amb;
              outer Real y1;
              outer Real x2;
              outer SI.Length d;
              outer SI.Length xLiner;
              outer SI.Length  xCFRP;
              outer SI.Length L;

        /****************** Tables *******************/
         Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            tableOnFile=true,
            table=[1,850,15,481; 2,2700,236,900; 3,1286,1.17,1578; 4,1374,1.14,1075],
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
        /****************** Initial equations *******************/
        initial equation
         T=T_amb;
        /****************** equations *******************/
        equation
        A=(d/2+xLiner+xCFRP+xLiner/10)*2*Modelica.Constants.pi*L+2*
        (d/2+xLiner+xCFRP)^2*Modelica.Constants.pi;

        Local_Radius = (d/2+xLiner+xCFRP+xLiner/10); // xLiner/10 accounts for the half thermal resistance in "InnerWallcell"

        tank_prop.u=y1;
        tank_prop.y[1]=rho;
        tank_prop.y[2]=k;
        tank_prop.y[3]=cp;

         R = dx / (A*k);

        portA.Q = (portA.T-T) / (R/2);
        //Q=h_o*(T_amb-T)*A;

         A*dx*rho*cp*der(T) = portA.Q + Q;

        Q_net = portA.Q + Q;

        if Adiabatic_Wall == true then
          Q = 0;

        else

        Q=h_o*(T_amb-T)*A;

        end if;

         annotation (Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},
                    {100,100}}),      graphics), Icon(coordinateSystem(
                 preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
               graphics={Rectangle(extent={{-90,60},{90,-60}}, lineColor={0,
                     0,255}), Text(
                 extent={{-78,22},{78,-14}},
                 lineColor={0,0,255},
                  textString="Outer_Wall")}));
        end OuterWallCell;

        model LinerCell
          import SI = Modelica.SIunits;

          /******************** Connectors *****************************/
         Ports.HeatFlow portA "Heat flow is pos when added to the wall"
            annotation (Placement(transformation(extent={{-10,50},{10,70}},
                  rotation=0)));
          Ports.HeatFlow portB "Heat flow is pos when subtracted from the wall"
            annotation (Placement(transformation(extent={{-10,-70},{10,-50}},
                  rotation=0)),HideResult=true);

        /****************** variables *******************/
          SI.Temperature T "Temperature of the wall element";
          SI.ThermalResistance R;
          SI.ThermalConductivity k;
          SI.SpecificHeatCapacity cp;
          SI.Density rho;
          SI.Area A "heat flow surface area, dz*dy";
          SI.Length dx=x1;
        SI.Length Radius_local;
        /****************** General parameters *******************/
         outer SI.Temperature T_amb;
         outer Real x1;
         outer Real y1;
         outer Real t1;
         outer SI.Length d;
         outer SI.Length L;

        /****************** Tables *******************/
        Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            table=[1,850,15,481; 2,2700,167,1106; 3,1286,1.17,1578; 4,1374,1.14,
                1075],
            tableOnFile=true,
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-74,46},{-54,66}})));

        /****************** Initial equations *******************/
        initial equation
          T=T_amb;
        /****************** equations *******************/
        equation
        A=((portA.Counter)*x1+d/2)*Modelica.Constants.pi*2*L+2*
        ((portA.Counter)*x1+d/2)^2*Modelica.Constants.pi;

        Radius_local = (portA.Counter)*x1+d/2;

        tank_prop.u=y1;
        tank_prop.y[1]=rho;
        tank_prop.y[2]=k;
        tank_prop.y[3]=cp;

          R = dx / (A*k);
          portA.Q = (portA.T-T) / (R/2);
          portB.Q = (portB.T-T) / (R/2);

          A*dx*rho*cp*der(T) = portA.Q + portB.Q;

        portA.Counter+1=portB.Counter;
        portA.P=portB.P;
          annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}),       graphics), Icon(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics={Rectangle(extent={{-90,60},{90,-60}}, lineColor={0,
                      0,255}), Text(
                  extent={{-82,36},{84,-38}},
                  lineColor={0,0,255},
                  textString="Liner")}));
        end LinerCell;

        model TankCell
          import SI = Modelica.SIunits;

          /******************** Connectors *****************************/
         Ports.HeatFlow portA "Heat flow is pos when added to the wall"
            annotation (Placement(transformation(extent={{-10,50},{10,70}},
                  rotation=0)));
         Ports.HeatFlow portB "Heat flow is pos when subtracted from the wall"
            annotation (Placement(transformation(extent={{-10,-70},{10,-50}},
                  rotation=0)),HideResult=true);

          /****************** variables *******************/
          SI.Temperature T "Temperature of the wall element";
          SI.ThermalResistance R;
          SI.ThermalConductivity k;
          SI.SpecificHeatCapacity cp;
          SI.Density rho;

          SI.Area A "heat flow surface area, dz*dy";
          SI.Length dx=x2;
          SI.Length Local_Radius;

          /****************** General parameters *******************/
          outer SI.Temperature T_amb;
          outer SI.Length x2;
          outer SI.Length x1;
          outer Real y2;
          outer SI.Length d;
          outer SI.Length L;
          outer Real t1;
          outer Real t2;

        /****************** Tables *******************/
        Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
            table=[1,850,15,481; 2,2700,236,900; 3,1286,1.17,1578; 4,1500,0.5,940],
            tableOnFile=true,
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-76,46},{-56,66}})));
          /****************** Initial equations *******************/
        initial equation
          T=T_amb;
          /****************** equations *******************/
        equation
        A=((portA.Counter)*x2+dx/2+x1*t1+d/2)*Modelica.Constants.pi*2*L+2*
        ((portA.Counter)*x2+dx/2+x1*t1+d/2)^2*Modelica.Constants.pi;

        Local_Radius = (portA.Counter*x2)+ dx/2 +x1*t1+d/2;

        tank_prop.u=y2;
        tank_prop.y[1]=rho;
        tank_prop.y[2]=k;
        tank_prop.y[3]=cp;

          R = dx / (A*k);
          portA.Q = (portA.T-T) / (R/2);
          portB.Q = (portB.T-T) / (R/2);

        A*dx*rho*cp*der(T) = portA.Q + portB.Q;
        portA.P=portB.P;
        if portA.Counter>=t2-0.5 then
          portB.Counter =0;
          else
        portA.Counter+1=portB.Counter;
        end if;
          annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}),       graphics), Icon(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics={Rectangle(extent={{-90,60},{90,-60}}, lineColor={0,
                      0,255}), Text(
                  extent={{-68,26},{66,-32}},
                  lineColor={0,0,255},
                  textString="CFRP")}));
        end TankCell;

        model TankCell_Inner
          import SI = Modelica.SIunits;

          /******************** Connectors *****************************/
         Ports.HeatFlow portA "Heat flow is pos when added to the wall"
            annotation (Placement(transformation(extent={{-10,50},{10,70}},
                  rotation=0)));
         Ports.HeatFlow portB "Heat flow is pos when subtracted from the wall"
            annotation (Placement(transformation(extent={{-10,-70},{10,-50}},
                  rotation=0)),HideResult=true);

          /****************** variables *******************/
          SI.Temperature T "Temperature of the wall element";
          SI.ThermalResistance R;
          SI.ThermalConductivity k;
          SI.SpecificHeatCapacity cp;
          SI.Density rho;

          SI.Area A "heat flow surface area, dz*dy";
          SI.Length dx=x2;
          SI.Length Local_Radius;

          /****************** General parameters *******************/
          outer SI.Temperature T_amb;
          outer SI.Length x2;
          outer SI.Length x1;
          outer Real y2;
          outer SI.Length d;
          outer SI.Length L;
          outer Real t1;
          outer Real t2;

        /****************** Tables *******************/
        Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
            table=[1,850,15,481; 2,2700,236,900; 3,1286,1.17,1578; 4,1500,0.5,940],
            tableOnFile=true,
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-76,46},{-56,66}})));
          /****************** Initial equations *******************/
        initial equation
          T=T_amb;

          /****************** equations *******************/
        equation
        A=(dx/2+x1*t1+d/2)*Modelica.Constants.pi*2*L+2*
        (dx/2+x1*t1+d/2)^2*Modelica.Constants.pi;

        Local_Radius = (dx/2+x1*t1+d/2);

        tank_prop.u=y2;
        tank_prop.y[1]=rho;
        tank_prop.y[2]=k;
        tank_prop.y[3]=cp;

          R = dx / (A*k);
          portA.Q = (portA.T-T) / (R/2);
          portB.Q = (portB.T-T) / (R/2);

        A*dx*rho*cp*der(T) = portA.Q + portB.Q;
        portA.P=portB.P;

        portA.Counter-portA.Counter +1=portB.Counter; // Restart portB.Counter = 1 for CFRP

          annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}),       graphics), Icon(coordinateSystem(
                  preserveAspectRatio=true,  extent={{-100,-100},{100,100}}),
                graphics={Rectangle(extent={{-90,60},{90,-60}}, lineColor={0,
                      0,255}), Text(
                  extent={{-68,26},{66,-32}},
                  lineColor={0,0,255},
                  textString="CFRP"),
                Text(
                  extent={{-88,54},{-16,20}},
                  lineColor={0,0,255},
                  textString="Inner")}));
        end TankCell_Inner;

        model TubeCell
          import SI = Modelica.SIunits;
          import C = Modelica.Constants;

        /******************** Connectors *****************************/
            Ports.TemperaturePort                     HT1
            annotation (Placement(transformation(extent={{-56,60},{-36,80}}),
                iconTransformation(extent={{-56,60},{-36,80}})));
          Ports.TemperaturePort                     HT2
            annotation (Placement(transformation(extent={{40,60},{60,80}})));

        /****************** General parameters *******************/
         outer SI.MassFlowRate m_flow;
         outer SI.SpecificHeatCapacity cp_h2;
         outer SI.CoefficientOfHeatTransfer h;
         outer SI.CoefficientOfHeatTransfer h_amb;
         outer parameter SI.Length d_i;
         outer parameter SI.Length d_o;
         outer parameter SI.Temperature T_amb;
         outer SI.Length L;
        parameter SI.Length dx=(d_o-d_i)/12;

        /****************** variables *******************/
        flow SI.HeatFlowRate[8] Q;
        SI.Temperature[7] T;
        Real R;
        SI.SpecificHeatCapacity cp;
        SI.Density rho;
        SI.Conductivity k;

        /****************** Tables *******************/
        Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            table=[1,850,15,481; 2,2700,167,1106; 3,1286,1.17,1578; 4,1374,1.14,
                1075],
            tableOnFile=true,
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-72,-10},{-52,10}})));
        /****************** Initial equations *******************/
        initial equation
          T[1]=T_amb;
          T[3]=T_amb;
          T[5]=T_amb;
          T[7]=T_amb;

        /****************** equations *******************/
        equation

          tank_prop.u=1;
          tank_prop.y[1]=rho;
          tank_prop.y[2]=k;
          tank_prop.y[3]=cp;

          R=dx/((d_i*C.pi*L)*k);
        //Between hydrogen and inner wall
          Q[1]=(HT1.T-T[1])/(1/(h*(d_i*C.pi*L)));
          Q[2]=(T[2]-T[1])/(R/2);
          (d_i*C.pi*L)*dx*rho*cp*der(T[1]) = Q[1] + Q[2];
        //Between two wall volumes
          -Q[2]=Q[3];
          Q[3]=(T[2]-T[3])/(R/2);
          Q[4]=(T[4]-T[3])/(R/2);
          ((d_i/2+dx*2)*2*C.pi*L)*2*dx*rho*cp*der(T[3]) = Q[3] + Q[4];
        //Between two wall volumes
          -Q[4]=Q[5];
          Q[5]=(T[4]-T[5])/(R/2);
          Q[6]=(T[6]-T[5])/(R/2);
          ((d_i/2+dx*4)*2*C.pi*L)*2*dx*rho*cp*der(T[5]) = Q[5] + Q[6];
        //Between outer wall volume and ambient
          -Q[6]=Q[7];
          Q[7]=(T[6]-T[7])/(R/2);
          Q[8]=h_amb*(T_amb-T[7])*((d_i/2+dx*6)*2*C.pi*L);
          ((d_i/2+dx*6)*2*C.pi*L)*dx*rho*cp*der(T[7]) = Q[7] + Q[8];

          -Q[1]=cp_h2*m_flow*(HT1.T-HT2.T);

          annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                    -100},{100,100}}), graphics), Icon(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics={Rectangle(extent={{-84,70},{96,-50}}, lineColor={0,
                      0,255}), Text(
                  extent={{-56,24},{66,-28}},
                  lineColor={0,0,255},
                  textString="Tube cell
")}));
        end TubeCell;

        model Liner5Pieces
        //Import of models
          LinerCell
                liner
            annotation (Placement(transformation(extent={{-10,50},{10,70}})));
          LinerCell
                liner1
            annotation (Placement(transformation(extent={{-10,24},{10,44}})));
          LinerCell
                liner2
            annotation (Placement(transformation(extent={{-10,0},{10,20}})));
          LinerCell
                liner3
            annotation (Placement(transformation(extent={{-10,-24},{10,-4}})));
          LinerCell
                liner4
            annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
          Ports.HeatFlow heatFlow
            annotation (Placement(transformation(extent={{-10,76},{10,96}})));
         Ports.HeatFlow heatFlow1
            annotation (Placement(transformation(extent={{-10,-80},{10,-60}})));
        equation
          connect(liner.portB, liner1.portA) annotation (Line(
              points={{0,54},{0,40}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(liner1.portB, liner2.portA) annotation (Line(
              points={{0,28},{0,16}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(liner2.portB, liner3.portA) annotation (Line(
              points={{0,4},{0,-8}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(liner3.portB, liner4.portA) annotation (Line(
              points={{0,-20},{0,-34}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(liner.portA, heatFlow) annotation (Line(
              points={{0,66},{0,86}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(liner4.portB, heatFlow1) annotation (Line(
              points={{0,-46},{0,-70}},
              color={0,0,255},
              smooth=Smooth.None));
          annotation (Diagram(graphics));
        end Liner5Pieces;

        model Tank10Pieces
        //import models
          TankCell
               cFRP1
            annotation (Placement(transformation(extent={{-50,12},{-30,32}})));
          TankCell
               cFRP2
            annotation (Placement(transformation(extent={{-50,-12},{-30,8}})));
          TankCell
               cFRP3
            annotation (Placement(transformation(extent={{-50,-36},{-30,-16}})));
          TankCell
               cFRP4
            annotation (Placement(transformation(extent={{-50,-60},{-30,-40}})));
          TankCell
               cFRP5
            annotation (Placement(transformation(extent={{-10,36},{10,56}})));
          TankCell
               cFRP6
            annotation (Placement(transformation(extent={{-10,12},{10,32}})));
          TankCell
               cFRP7
            annotation (Placement(transformation(extent={{-10,-12},{10,8}})));
          TankCell
               cFRP8
            annotation (Placement(transformation(extent={{-10,-36},{10,-16}})));
          TankCell
               cFRP9
            annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));
          Ports.HeatFlow heatFlow
            annotation (Placement(transformation(extent={{-10,74},{10,94}})));
          Ports.HeatFlow heatFlow1
            annotation (Placement(transformation(extent={{-10,-88},{10,-68}})));
          TankCell_Inner tankCell_Inner
            annotation (Placement(transformation(extent={{-50,36},{-30,56}})));
        equation
          connect(cFRP1.portB, cFRP2.portA) annotation (Line(
              points={{-40,16},{-40,4}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP2.portB, cFRP3.portA) annotation (Line(
              points={{-40,-8},{-40,-20}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP3.portB, cFRP4.portA) annotation (Line(
              points={{-40,-32},{-40,-44}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP4.portB, cFRP5.portA) annotation (Line(
              points={{-40,-56},{-40,-64},{-18,-64},{-18,60},{0,60},{0,52}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP5.portB, cFRP6.portA) annotation (Line(
              points={{0,40},{0,28}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP6.portB, cFRP7.portA) annotation (Line(
              points={{0,16},{0,4}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP7.portB, cFRP8.portA) annotation (Line(
              points={{0,-8},{0,-20}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP8.portB, cFRP9.portA) annotation (Line(
              points={{0,-32},{0,-44}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP9.portB, heatFlow1) annotation (Line(
              points={{0,-56},{0,-56},{0,-78}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(tankCell_Inner.portB, cFRP1.portA) annotation (Line(
              points={{-40,40},{-40,28}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tankCell_Inner.portA, heatFlow) annotation (Line(
              points={{-40,52},{-40,84},{0,84}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          annotation (Diagram(graphics));
        end Tank10Pieces;

        model Tube5Pieces

          TubeCell   tubeHeatTransfer3Testing
            annotation (Placement(transformation(extent={{-60,4},{-40,24}})));
          TubeCell   tubeHeatTransfer3Testing1
            annotation (Placement(transformation(extent={{-34,4},{-14,24}})));
          TubeCell   tubeHeatTransfer3Testing2
            annotation (Placement(transformation(extent={{-8,4},{12,24}})));
          TubeCell  tubeHeatTransfer3Testing3
            annotation (Placement(transformation(extent={{18,4},{38,24}})));
          TubeCell  tubeHeatTransfer3Testing4
            annotation (Placement(transformation(extent={{44,4},{64,24}})));

          Ports.TemperaturePort                     HT1
            annotation (Placement(transformation(extent={{-72,50},{-52,70}}),
                iconTransformation(extent={{-72,50},{-52,70}})));

          Ports.TemperaturePort                     HT2
            annotation (Placement(transformation(extent={{48,50},{68,70}}),
                iconTransformation(extent={{48,50},{68,70}})));
        equation
          connect(HT1,tubeHeatTransfer3Testing. HT1) annotation (Line(
              points={{-62,60},{-54.6,60},{-54.6,21}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tubeHeatTransfer3Testing1.HT2,tubeHeatTransfer3Testing2. HT1)
            annotation (Line(
              points={{-19,21},{-10,22},{-2.6,21}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tubeHeatTransfer3Testing2.HT2,tubeHeatTransfer3Testing3. HT1)
            annotation (Line(
              points={{7,21},{7,18.5},{23.4,18.5},{23.4,21}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tubeHeatTransfer3Testing3.HT2,tubeHeatTransfer3Testing4. HT1)
            annotation (Line(
              points={{33,21},{49.4,21}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tubeHeatTransfer3Testing4.HT2, HT2) annotation (Line(
              points={{59,21},{59,28.5},{58,28.5},{58,60}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tubeHeatTransfer3Testing1.HT1, tubeHeatTransfer3Testing.HT2)
            annotation (Line(
              points={{-28.6,21},{-45,21}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                    -100},{100,100}}), graphics), Icon(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics));
        end Tube5Pieces;
      end WallPieces_Area_middle;

      package WallPieces_PCM_Area_middle

        model InnerWallCell

          import SI = Modelica.SIunits;
        /*********************** Thermodynamic property call ***********************************/

                 /*    replaceable package Medium = CoolProp2Modelica.Media.Hydrogen (
        onePhase=true) constrainedby Modelica.Media.Interfaces.PartialMedium 
    annotation (choicesAllMatching=true);
*/

         replaceable package Medium =
            ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
            Modelica.Media.Interfaces.PartialMedium;
            Medium.ThermodynamicState medium;

        /******************** Connectors *****************************/
          Ports.HeatFlow2 portA "Heat flow is pos when added to the wall"
            annotation (Placement(transformation(extent={{-10,50},{10,70}},
                  rotation=0)));
          Ports.HeatFlow portB "Heat flow is pos when subtracted from the wall"
            annotation (Placement(transformation(extent={{-10,-70},{10,-50}},
                  rotation=0)),HideResult=true);

        /****************** variables *******************/
          SI.Temperature T "Temperature of the wall element";
          SI.Pressure p "Pressure of the hydrogen";
          SI.ThermalResistance R "Thermal resistance of liner";
          SI.CoefficientOfHeatTransfer h "Heat transfer coefficient";
          SI.ThermalConductivity k "Thermal conductivity of liner";
          SI.SpecificHeatCapacity cp "Specific heat capacity of liner";
          SI.Density rho "Density of liner";
          Real Tau "Dimenasionless time";
          Real Ra "Dimensionless hydrogen properties number";
          Real beta "Thermal expansion coefficient";
          Real v "kinematic viscosity";
          Real a "thermal diffusivity of hydrogen";
          Real Nu "Dimensionless heat transfer number";

        Real h_TEST;

        /****************** General parameters *******************/
          constant Real g=9.82 "acceleration due to gravity";
          SI.Length dx=x1/2;

              outer SI.Length d "inside diameter of cylinder";
              outer SI.Temperature T_amb "Ambient temperature";
              outer Real y1;
              outer SI.CoefficientOfHeatTransfer  h_i;
              outer SI.Area A;
              outer SI.Length x1;

        /****************** Tables *******************/
          Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            tableOnFile=true,
            table=[1,850,15,481; 2,2700,236,900; 3,1286,1.17,1578; 4,1374,1.14,1075],
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-80,0},{-60,20}})));

        /****************** Initial equations *******************/
        initial equation
          T=T_amb;

        /****************** equations *******************/
        equation
        medium=Medium.setState_pT(p, T);
        tank_prop.u=y1;
        tank_prop.y[1]=rho;
        tank_prop.y[2]=k;
        tank_prop.y[3]=cp;

         Medium.thermalConductivity(medium)=a;
         Medium.dynamicViscosity(medium)/Medium.density(medium)=v;

        //calculation of rayleighs number taken from 'Natural convection cooling of
        //rectangular and cylindrical containers' by Wenxian Lin, S.W. Armfield
        Ra=g*beta*d^3*Medium.specificHeatCapacityCp(medium)*
        Medium.density(medium)^2*abs(portA.T-T)/(v*k);
        beta=1/portA.T;
        Tau=time/(d^2/(a*Ra^(1/2)));
        Nu=0.104*Ra^(0.352);

        h_TEST = Nu*k/d;

        /*
if portA.m_flow < 0.00 then
  h=Nu*k/d;
elseif portA.m_flow > 0.00 then
  h=h_i;
else
  h=50;
end if;
*/

        if portA.m_flow > 0 then
        h = 150; // equal to "h_charging" in "HeatTransferTank"
        else

          h = 50;

        end if;

        R = dx / (A*k);

          portA.Q = (portA.T-T) / (1/(h*A));
          portB.Q = (portB.T-T) / R; // instead of R/2, because there is only one conduction resistance and dx = x1/2 already

          A*dx*rho*cp*der(T) = portA.Q + portB.Q;

        portA.P=p;
        portA.P=portB.P;
        portA.Counter+1=portB.Counter;
          annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}),       graphics), Icon(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics={Rectangle(extent={{-90,60},{90,-60}}, lineColor={0,
                      0,255}), Text(
                  extent={{-86,16},{86,-4}},
                  lineColor={0,0,255},
                  textString="Inner_wall_discharging")}));
        end InnerWallCell;

        model OuterWallCell "Outer wall with natural convection given by h_o"
          import SI = Modelica.SIunits;

        /******************** Connectors *****************************/
         Ports.HeatFlow portA "Heat flow is pos when added to the wall"
           annotation (Placement(transformation(extent={{-10,50},{10,70}},
                 rotation=0)));

        /****************** variables *******************/
         SI.Temperature T "Temperature of the wall element";
         SI.ThermalResistance R;
         SI.ThermalConductivity k;
         SI.SpecificHeatCapacity cp;
         SI.Density rho;
         SI.Area  A "heat flow surface area, dz*dy";
         SI.Length dx=x2/2;
         SI.HeatFlowRate Q;
          SI.HeatFlowRate Q_net;
        SI.Length Radius_Local;
        outer parameter Boolean Adiabatic_Wall;
        /****************** General parameters *******************/
              outer parameter SI.CoefficientOfHeatTransfer h_o;
              outer SI.Temperature T_amb;
              outer Real y1;
              outer Real x2;
              outer SI.Length d;
              outer SI.Length xLiner;
              outer SI.Length  xCFRP;
              outer SI.Length L;
         outer SI.Length t_PCM;
         outer SI.Length dx_PCM;

        /****************** Tables *******************/
         Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            tableOnFile=true,
            table=[1,850,15,481; 2,2700,236,900; 3,1286,1.17,1578; 4,1374,1.14,1075],
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
        /****************** Initial equations *******************/
        initial equation
         T=T_amb;
        /****************** equations *******************/
        equation
        A=(d/2 + t_PCM + dx_PCM/2 + xLiner  + xCFRP)*2*Modelica.Constants.pi*L+2*
        (d/2 + t_PCM + dx_PCM/2 + xLiner  + xCFRP)^2*Modelica.Constants.pi;

        Radius_Local = d/2 + t_PCM + dx_PCM/2 + xLiner  + xCFRP;

        tank_prop.u=y1;
        tank_prop.y[1]=rho;
        tank_prop.y[2]=k;
        tank_prop.y[3]=cp;

         R = dx / (A*k);

        portA.Q = (portA.T-T) / (R/2);
        //Q=h_o*(T_amb-T)*A;

         A*dx*rho*cp*der(T) = portA.Q + Q;

        Q_net = portA.Q + Q;

        if Adiabatic_Wall == true then
          Q = 0;

        else

        Q=h_o*(T_amb-T)*A;

        end if;

         annotation (Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},
                    {100,100}}),      graphics), Icon(coordinateSystem(
                 preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
               graphics={Rectangle(extent={{-90,60},{90,-60}}, lineColor={0,
                     0,255}), Text(
                 extent={{-78,22},{78,-14}},
                 lineColor={0,0,255},
                  textString="Outer_Wall")}));
        end OuterWallCell;

        model LinerCell
          import SI = Modelica.SIunits;

          /******************** Connectors *****************************/
         Ports.HeatFlow portA "Heat flow is pos when added to the wall"
            annotation (Placement(transformation(extent={{-10,50},{10,70}},
                  rotation=0)));
          Ports.HeatFlow portB "Heat flow is pos when subtracted from the wall"
            annotation (Placement(transformation(extent={{-10,-70},{10,-50}},
                  rotation=0)),HideResult=true);

        /****************** variables *******************/
          SI.Temperature T "Temperature of the wall element";
          SI.ThermalResistance R;
          SI.ThermalConductivity k;
          SI.SpecificHeatCapacity cp;
          SI.Density rho;
          SI.Area A "heat flow surface area, dz*dy";
          SI.Length dx=x1;
        SI.Length Radius_local;
        /****************** General parameters *******************/
         outer SI.Temperature T_amb;
         outer Real x1;
         outer Real y1;
         outer Real t1;
         outer SI.Length d;
         outer SI.Length L;
         outer SI.Length t_PCM;
         outer SI.Length dx_PCM;

        /****************** Tables *******************/
        Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            table=[1,850,15,481; 2,2700,167,1106; 3,1286,1.17,1578; 4,1374,1.14,
                1075],
            tableOnFile=true,
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-74,46},{-54,66}})));

        /****************** Initial equations *******************/
        initial equation
          T=T_amb;
        /****************** equations *******************/
        equation
        A=(d/2 + t_PCM + dx_PCM/2+ x1/2 + portA.Counter*x1)*Modelica.Constants.pi*2*L+2*
        (d/2 + t_PCM + dx_PCM/2+ x1/2 + portA.Counter*x1)^2*Modelica.Constants.pi;

        Radius_local = (d/2 + t_PCM + dx_PCM/2+ x1/2 + portA.Counter*x1);

        tank_prop.u=y1;
        tank_prop.y[1]=rho;
        tank_prop.y[2]=k;
        tank_prop.y[3]=cp;

          R = dx / (A*k);
          portA.Q = (portA.T-T) / (R/2);
          portB.Q = (portB.T-T) / (R/2);

          A*dx*rho*cp*der(T) = portA.Q + portB.Q;

        portA.Counter+1=portB.Counter;
        portA.P=portB.P;
          annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}),       graphics), Icon(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics={Rectangle(extent={{-90,60},{90,-60}}, lineColor={0,
                      0,255}), Text(
                  extent={{-82,36},{84,-38}},
                  lineColor={0,0,255},
                  textString="Liner")}));
        end LinerCell;

        model LinerCell_Inner
          import SI = Modelica.SIunits;

          /******************** Connectors *****************************/
         Ports.HeatFlow portA "Heat flow is pos when added to the wall"
            annotation (Placement(transformation(extent={{-10,50},{10,70}},
                  rotation=0)));
          Ports.HeatFlow portB "Heat flow is pos when subtracted from the wall"
            annotation (Placement(transformation(extent={{-10,-70},{10,-50}},
                  rotation=0)),HideResult=true);

        /****************** variables *******************/
          SI.Temperature T "Temperature of the wall element";
          SI.ThermalResistance R;
          SI.ThermalConductivity k;
          SI.SpecificHeatCapacity cp;
          SI.Density rho;
          SI.Area A "heat flow surface area, dz*dy";
          SI.Length dx=x1;
        SI.Length Radius_local;
        /****************** General parameters *******************/
         outer SI.Temperature T_amb;
         outer Real x1;
         outer Real y1;
         outer Real t1;
         outer SI.Length d;
         outer SI.Length L;
         outer SI.Length t_PCM;
         outer SI.Length dx_PCM;
        /****************** Tables *******************/
        Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            table=[1,850,15,481; 2,2700,167,1106; 3,1286,1.17,1578; 4,1374,1.14,
                1075],
            tableOnFile=true,
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-74,46},{-54,66}})));

        /****************** Initial equations *******************/
        initial equation
          T=T_amb;
        /****************** equations *******************/
        equation
        A=(d/2 + t_PCM + dx_PCM/2+ x1/2)*Modelica.Constants.pi*2*L+2*
        (d/2 + t_PCM + dx_PCM/2+ x1/2)^2*Modelica.Constants.pi;

        // A=(dx/2+x1*t1+d/2)*Modelica.Constants.pi*2*L+2*(dx/2+x1*t1+d/2)^2*Modelica.Constants.pi;
        Radius_local = d/2 + t_PCM + dx_PCM/2 + x1/2;

        tank_prop.u=y1;
        tank_prop.y[1]=rho;
        tank_prop.y[2]=k;
        tank_prop.y[3]=cp;

          R = dx / (A*k);
          portA.Q = (portA.T-T) / (R/2);
          portB.Q = (portB.T-T) / (R/2);

          A*dx*rho*cp*der(T) = portA.Q + portB.Q;

        portA.Counter - portA.Counter +1=portB.Counter; // reset portB.Counter = 1
        portA.P=portB.P;
          annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}),       graphics), Icon(coordinateSystem(
                  preserveAspectRatio=true,  extent={{-100,-100},{100,100}}),
                graphics={Rectangle(extent={{-90,60},{90,-60}}, lineColor={0,
                      0,255}), Text(
                  extent={{-82,36},{84,-38}},
                  lineColor={0,0,255},
                  textString="Liner"),
                Text(
                  extent={{-100,62},{-2,24}},
                  lineColor={0,0,255},
                  textString="Inner")}));
        end LinerCell_Inner;

        model TankCell
          import SI = Modelica.SIunits;

          /******************** Connectors *****************************/
         Ports.HeatFlow portA "Heat flow is pos when added to the wall"
            annotation (Placement(transformation(extent={{-10,50},{10,70}},
                  rotation=0)));
         Ports.HeatFlow portB "Heat flow is pos when subtracted from the wall"
            annotation (Placement(transformation(extent={{-10,-70},{10,-50}},
                  rotation=0)),HideResult=true);

          /****************** variables *******************/
          SI.Temperature T "Temperature of the wall element";
          SI.ThermalResistance R;
          SI.ThermalConductivity k;
          SI.SpecificHeatCapacity cp;
          SI.Density rho;

          SI.Area A "heat flow surface area, dz*dy";
          SI.Length dx=x2;
          SI.Length Radius_Local;

          /****************** General parameters *******************/
          outer SI.Temperature T_amb;
          outer SI.Length x2;
          outer SI.Length x1;
          outer Real y2;
          outer SI.Length d;
          outer SI.Length L;
          outer Real t1;
          outer Real t2;
         outer SI.Length t_PCM;
         outer SI.Length dx_PCM;

        /****************** Tables *******************/
        Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
            table=[1,850,15,481; 2,2700,236,900; 3,1286,1.17,1578; 4,1500,0.5,940],
            tableOnFile=true,
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-76,46},{-56,66}})));
          /****************** Initial equations *******************/
        initial equation
          T=T_amb;
          /****************** equations *******************/
        equation
        A=(d/2 + t_PCM + dx_PCM/2 + (t1-0.5)*x1 + x2/2 + portA.Counter*x2)*Modelica.Constants.pi*2*L+2*
        (d/2 + t_PCM + dx_PCM/2 + (t1-0.5)*x1 + x2/2 + portA.Counter*x2)^2*Modelica.Constants.pi;

        Radius_Local = (d/2 + t_PCM + dx_PCM/2 + (t1-0.5)*x1 + x2/2 + portA.Counter*x2);

        tank_prop.u=y2;
        tank_prop.y[1]=rho;
        tank_prop.y[2]=k;
        tank_prop.y[3]=cp;

          R = dx / (A*k);
          portA.Q = (portA.T-T) / (R/2);
          portB.Q = (portB.T-T) / (R/2);

        A*dx*rho*cp*der(T) = portA.Q + portB.Q;
        portA.P=portB.P;
        if portA.Counter>=t2-0.5 then
          portB.Counter =0;
          else
        portA.Counter+1=portB.Counter;
        end if;
          annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}),       graphics), Icon(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics={Rectangle(extent={{-90,60},{90,-60}}, lineColor={0,
                      0,255}), Text(
                  extent={{-68,26},{66,-32}},
                  lineColor={0,0,255},
                  textString="CFRP")}));
        end TankCell;

        model TankCell_Inner
          import SI = Modelica.SIunits;

          /******************** Connectors *****************************/
         Ports.HeatFlow portA "Heat flow is pos when added to the wall"
            annotation (Placement(transformation(extent={{-10,50},{10,70}},
                  rotation=0)));
         Ports.HeatFlow portB "Heat flow is pos when subtracted from the wall"
            annotation (Placement(transformation(extent={{-10,-70},{10,-50}},
                  rotation=0)),HideResult=true);

          /****************** variables *******************/
          SI.Temperature T "Temperature of the wall element";
          SI.ThermalResistance R;
          SI.ThermalConductivity k;
          SI.SpecificHeatCapacity cp;
          SI.Density rho;

          SI.Area A "heat flow surface area, dz*dy";
          SI.Length dx=x2;
          SI.Length Radius_Local;

          /****************** General parameters *******************/
          outer SI.Temperature T_amb;
          outer SI.Length x2;
          outer SI.Length x1;
          outer Real y2;
          outer SI.Length d;
          outer SI.Length L;
          outer Real t1;
          outer Real t2;
          outer SI.Length t_PCM;
          outer SI.Length dx_PCM;

        /****************** Tables *******************/
        Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
            table=[1,850,15,481; 2,2700,236,900; 3,1286,1.17,1578; 4,1500,0.5,940],
            tableOnFile=true,
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-76,46},{-56,66}})));
          /****************** Initial equations *******************/
        initial equation
          T=T_amb;

          /****************** equations *******************/
        equation
        A=(d/2 + t_PCM + dx_PCM/2 + (t1-0.5)*x1 + x2/2)*Modelica.Constants.pi*2*L+2*
        (d/2 + t_PCM + dx_PCM/2 + (t1-0.5)*x1 + x2/2)^2*Modelica.Constants.pi;

        Radius_Local = (d/2 + t_PCM + dx_PCM/2 + (t1-0.5)*x1 + x2/2);

        tank_prop.u=y2;
        tank_prop.y[1]=rho;
        tank_prop.y[2]=k;
        tank_prop.y[3]=cp;

          R = dx / (A*k);
          portA.Q = (portA.T-T) / (R/2);
          portB.Q = (portB.T-T) / (R/2);

        A*dx*rho*cp*der(T) = portA.Q + portB.Q;
        portA.P=portB.P;

        portA.Counter-portA.Counter +1=portB.Counter; // Restart portB.Counter = 1 for CFRP

          annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}),       graphics), Icon(coordinateSystem(
                  preserveAspectRatio=true,  extent={{-100,-100},{100,100}}),
                graphics={Rectangle(extent={{-90,60},{90,-60}}, lineColor={0,
                      0,255}), Text(
                  extent={{-68,26},{66,-32}},
                  lineColor={0,0,255},
                  textString="CFRP"),
                Text(
                  extent={{-88,54},{-16,20}},
                  lineColor={0,0,255},
                  textString="Inner")}));
        end TankCell_Inner;

        model TubeCell
          import SI = Modelica.SIunits;
          import C = Modelica.Constants;

        /******************** Connectors *****************************/
            Ports.TemperaturePort                     HT1
            annotation (Placement(transformation(extent={{-56,60},{-36,80}}),
                iconTransformation(extent={{-56,60},{-36,80}})));
          Ports.TemperaturePort                     HT2
            annotation (Placement(transformation(extent={{40,60},{60,80}})));

        /****************** General parameters *******************/
         outer SI.MassFlowRate m_flow;
         outer SI.SpecificHeatCapacity cp_h2;
         outer SI.CoefficientOfHeatTransfer h;
         outer SI.CoefficientOfHeatTransfer h_amb;
         outer parameter SI.Length d_i;
         outer parameter SI.Length d_o;
         outer parameter SI.Temperature T_amb;
         outer SI.Length L;
        parameter SI.Length dx=(d_o-d_i)/12;

        /****************** variables *******************/
        flow SI.HeatFlowRate[8] Q;
        SI.Temperature[7] T;
        Real R;
        SI.SpecificHeatCapacity cp;
        SI.Density rho;
        SI.Conductivity k;

        /****************** Tables *******************/
        Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
            tableName="properties",
            table=[1,850,15,481; 2,2700,167,1106; 3,1286,1.17,1578; 4,1374,1.14,
                1075],
            tableOnFile=true,
            fileName=
                "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
            annotation (Placement(transformation(extent={{-72,-10},{-52,10}})));
        /****************** Initial equations *******************/
        initial equation
          T[1]=T_amb;
          T[3]=T_amb;
          T[5]=T_amb;
          T[7]=T_amb;

        /****************** equations *******************/
        equation

          tank_prop.u=1;
          tank_prop.y[1]=rho;
          tank_prop.y[2]=k;
          tank_prop.y[3]=cp;

          R=dx/((d_i*C.pi*L)*k);
        //Between hydrogen and inner wall
          Q[1]=(HT1.T-T[1])/(1/(h*(d_i*C.pi*L)));
          Q[2]=(T[2]-T[1])/(R/2);
          (d_i*C.pi*L)*dx*rho*cp*der(T[1]) = Q[1] + Q[2];
        //Between two wall volumes
          -Q[2]=Q[3];
          Q[3]=(T[2]-T[3])/(R/2);
          Q[4]=(T[4]-T[3])/(R/2);
          ((d_i/2+dx*2)*2*C.pi*L)*2*dx*rho*cp*der(T[3]) = Q[3] + Q[4];
        //Between two wall volumes
          -Q[4]=Q[5];
          Q[5]=(T[4]-T[5])/(R/2);
          Q[6]=(T[6]-T[5])/(R/2);
          ((d_i/2+dx*4)*2*C.pi*L)*2*dx*rho*cp*der(T[5]) = Q[5] + Q[6];
        //Between outer wall volume and ambient
          -Q[6]=Q[7];
          Q[7]=(T[6]-T[7])/(R/2);
          Q[8]=h_amb*(T_amb-T[7])*((d_i/2+dx*6)*2*C.pi*L);
          ((d_i/2+dx*6)*2*C.pi*L)*dx*rho*cp*der(T[7]) = Q[7] + Q[8];

          -Q[1]=cp_h2*m_flow*(HT1.T-HT2.T);

          annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                    -100},{100,100}}), graphics), Icon(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics={Rectangle(extent={{-84,70},{96,-50}}, lineColor={0,
                      0,255}), Text(
                  extent={{-56,24},{66,-28}},
                  lineColor={0,0,255},
                  textString="Tube cell
")}));
        end TubeCell;

        model Liner5Pieces
        //Import of models
          LinerCell
                liner1
            annotation (Placement(transformation(extent={{-10,24},{10,44}})));
          LinerCell
                liner2
            annotation (Placement(transformation(extent={{-10,0},{10,20}})));
          LinerCell
                liner3
            annotation (Placement(transformation(extent={{-10,-24},{10,-4}})));
          LinerCell
                liner4
            annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
          Ports.HeatFlow heatFlow
            annotation (Placement(transformation(extent={{-10,76},{10,96}})));
         Ports.HeatFlow heatFlow1
            annotation (Placement(transformation(extent={{-10,-80},{10,-60}})));
          LinerCell_Inner linerCell_Inner
            annotation (Placement(transformation(extent={{-10,50},{10,70}})));
        equation
          connect(liner1.portB, liner2.portA) annotation (Line(
              points={{0,28},{0,16}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(liner2.portB, liner3.portA) annotation (Line(
              points={{0,4},{0,-8}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(liner3.portB, liner4.portA) annotation (Line(
              points={{0,-20},{0,-34}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(liner4.portB, heatFlow1) annotation (Line(
              points={{0,-46},{0,-70}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(heatFlow, linerCell_Inner.portA) annotation (Line(
              points={{0,86},{0,66}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(linerCell_Inner.portB, liner1.portA) annotation (Line(
              points={{0,54},{0,40}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          annotation (Diagram(graphics));
        end Liner5Pieces;

        model Tank10Pieces
        //import models
          TankCell
               cFRP1
            annotation (Placement(transformation(extent={{-50,12},{-30,32}})));
          TankCell
               cFRP2
            annotation (Placement(transformation(extent={{-50,-12},{-30,8}})));
          TankCell
               cFRP3
            annotation (Placement(transformation(extent={{-50,-36},{-30,-16}})));
          TankCell
               cFRP4
            annotation (Placement(transformation(extent={{-50,-60},{-30,-40}})));
          TankCell
               cFRP5
            annotation (Placement(transformation(extent={{-10,36},{10,56}})));
          TankCell
               cFRP6
            annotation (Placement(transformation(extent={{-10,12},{10,32}})));
          TankCell
               cFRP7
            annotation (Placement(transformation(extent={{-10,-12},{10,8}})));
          TankCell
               cFRP8
            annotation (Placement(transformation(extent={{-10,-36},{10,-16}})));
          TankCell
               cFRP9
            annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));
          Ports.HeatFlow heatFlow
            annotation (Placement(transformation(extent={{-10,74},{10,94}})));
          Ports.HeatFlow heatFlow1
            annotation (Placement(transformation(extent={{-10,-88},{10,-68}})));
          TankCell_Inner tankCell_Inner
            annotation (Placement(transformation(extent={{-50,36},{-30,56}})));
        equation
          connect(cFRP1.portB, cFRP2.portA) annotation (Line(
              points={{-40,16},{-40,4}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP2.portB, cFRP3.portA) annotation (Line(
              points={{-40,-8},{-40,-20}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP3.portB, cFRP4.portA) annotation (Line(
              points={{-40,-32},{-40,-44}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP4.portB, cFRP5.portA) annotation (Line(
              points={{-40,-56},{-40,-64},{-18,-64},{-18,60},{0,60},{0,52}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP5.portB, cFRP6.portA) annotation (Line(
              points={{0,40},{0,28}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP6.portB, cFRP7.portA) annotation (Line(
              points={{0,16},{0,4}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP7.portB, cFRP8.portA) annotation (Line(
              points={{0,-8},{0,-20}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP8.portB, cFRP9.portA) annotation (Line(
              points={{0,-32},{0,-44}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(cFRP9.portB, heatFlow1) annotation (Line(
              points={{0,-56},{0,-56},{0,-78}},
              color={0,0,255},
              smooth=Smooth.None));
          connect(tankCell_Inner.portB, cFRP1.portA) annotation (Line(
              points={{-40,40},{-40,28}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tankCell_Inner.portA, heatFlow) annotation (Line(
              points={{-40,52},{-40,84},{0,84}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          annotation (Diagram(graphics));
        end Tank10Pieces;

        model Tube5Pieces

          TubeCell   tubeHeatTransfer3Testing
            annotation (Placement(transformation(extent={{-60,4},{-40,24}})));
          TubeCell   tubeHeatTransfer3Testing1
            annotation (Placement(transformation(extent={{-34,4},{-14,24}})));
          TubeCell   tubeHeatTransfer3Testing2
            annotation (Placement(transformation(extent={{-8,4},{12,24}})));
          TubeCell  tubeHeatTransfer3Testing3
            annotation (Placement(transformation(extent={{18,4},{38,24}})));
          TubeCell  tubeHeatTransfer3Testing4
            annotation (Placement(transformation(extent={{44,4},{64,24}})));

          Ports.TemperaturePort                     HT1
            annotation (Placement(transformation(extent={{-72,50},{-52,70}}),
                iconTransformation(extent={{-72,50},{-52,70}})));

          Ports.TemperaturePort                     HT2
            annotation (Placement(transformation(extent={{48,50},{68,70}}),
                iconTransformation(extent={{48,50},{68,70}})));
        equation
          connect(HT1,tubeHeatTransfer3Testing. HT1) annotation (Line(
              points={{-62,60},{-54.6,60},{-54.6,21}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tubeHeatTransfer3Testing1.HT2,tubeHeatTransfer3Testing2. HT1)
            annotation (Line(
              points={{-19,21},{-10,22},{-2.6,21}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tubeHeatTransfer3Testing2.HT2,tubeHeatTransfer3Testing3. HT1)
            annotation (Line(
              points={{7,21},{7,18.5},{23.4,18.5},{23.4,21}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tubeHeatTransfer3Testing3.HT2,tubeHeatTransfer3Testing4. HT1)
            annotation (Line(
              points={{33,21},{49.4,21}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tubeHeatTransfer3Testing4.HT2, HT2) annotation (Line(
              points={{59,21},{59,28.5},{58,28.5},{58,60}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          connect(tubeHeatTransfer3Testing1.HT1, tubeHeatTransfer3Testing.HT2)
            annotation (Line(
              points={{-28.6,21},{-45,21}},
              color={0,0,0},
              thickness=0.5,
              smooth=Smooth.None));
          annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                    -100},{100,100}}), graphics), Icon(coordinateSystem(
                  preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics));
        end Tube5Pieces;
      end WallPieces_PCM_Area_middle;

      model HeatTransferTank_no_PCM

        import SI = Modelica.SIunits;

      /******************** Connectors *****************************/
        Ports.HeatFlow2 heatFlow
          annotation (Placement(transformation(extent={{-26,80},{-6,100}})));

      /****************** General parameters *******************/
      parameter Integer tank=4 "Tank type" annotation (Dialog(group= " Heat Transfer"),choices(
      choice=1 "Tyep 1 Steel",
      choice=2 "Tyep 2 Aluminium",
      choice=3 "Tyep 3 Aluminium liner",
      choice=4 "Type 4 Plastic liner"));

        inner parameter Boolean Charging = true
          "If true, tank is charging with given heat transfer coefficient" annotation(Dialog(group= " Heat Transfer"));
        parameter SI.CoefficientOfHeatTransfer  h_charging=150 "Charging, Heat transfer 
  coefficient inside HSS tank 150 w/m2K - 500w/m2K  
    according to monde (0 is adiabatic fuelling)"                                                             annotation(Dialog(group= " Heat Transfer"));

        parameter SI.CoefficientOfHeatTransfer  h_discharging=-1 "Discharging heat transfer 
  coefficient, if <0 then Daney relation
     is used, if >0 then the given number is used "                                                           annotation(Dialog(group= " Heat Transfer"));
        inner parameter SI.CoefficientOfHeatTransfer  h_o=8
          "Heat transfer coefficient outside the tanks (typical natural convection)"
                                                                                     annotation(Dialog(group= " Heat Transfer"));
       // inner parameter Real  T_amb=273;

       inner parameter SI.Length xLiner=0.003 "Thickness of liner" annotation(Dialog(group= "Geometry"));
       inner parameter SI.Length  xCFRP=0.022 "Thickness of wrapping/tank" annotation(Dialog(group= "Geometry"));

      protected
       inner parameter SI.Length x1=xLiner/(t1-0.5);
       inner parameter SI.Length x2=xCFRP/(t2-0.5);
       inner constant Real t1=5.5;
       inner constant Real t2=10.5;
       inner Real y1;
       inner Real y2;

      public
       inner SI.CoefficientOfHeatTransfer  h_i = h_charging;
       //(if Charging == true then h_charging else h_discharging);
      /*
 inner SI.Length   d = (if Area == false then AInner/(Modelica.Constants.pi)
 else dInner);
 inner SI.Length   L= (if Area == false then 0.8 else LInner);
 inner SI.Area   A = (if Area == false then AInner else 2*Modelica.Constants.pi*d/2*L+2*(d/2)^2*Modelica.Constants.pi);
 */
      inner parameter Boolean Adiabatic_Wall = true;

      /****************** equations *******************/

       WallPieces_Area_middle.InnerWallCell
                               wallCell_discharging
          annotation (Placement(transformation(extent={{-32,50},{-12,70}})));
        WallPieces_Area_middle.OuterWallCell
                               outer_wall
          annotation (Placement(transformation(extent={{-10,-38},{10,-18}})));
        WallPieces_Area_middle.Liner5Pieces wall_liner
          annotation (Placement(transformation(extent={{-4,14},{16,34}})));
        WallPieces_Area_middle.Tank10Pieces wall_CFRP
          annotation (Placement(transformation(extent={{-40,-16},{-20,4}})));

      equation
      //Deciding which calls to make in lookup tables in wallpieces for tank properties
       if tank==1 then
         y1=1;
         y2=1;
       elseif tank==2 then
         y1=2;
         y2=2;
       elseif tank==3 then
         y1=2;
         y2=4;
       elseif tank==4 then
         y1=3;
         y2=4;
       else
         y1=3;
         y2=4;
         end if;

        connect(heatFlow, wallCell_discharging.portA) annotation (Line(
            points={{-16,90},{-20,90},{-20,66},{-22,66}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(wallCell_discharging.portB, wall_liner.heatFlow) annotation (Line(
            points={{-22,54},{-22,38},{6,38},{6,32.6}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(wall_liner.heatFlow1, wall_CFRP.heatFlow) annotation (Line(
            points={{6,17},{6,8},{-30,8},{-30,2.4}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(wall_CFRP.heatFlow1, outer_wall.portA) annotation (Line(
            points={{-30,-13.8},{-14,-13.8},{-14,-16},{0,-16},{0,-22}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
                  -100},{100,100}}),
                            graphics), Icon(coordinateSystem(preserveAspectRatio=
                  false, extent={{-100,-100},{100,100}}),
                                            graphics={Bitmap(
                extent={{-90,-80},{96,74}},
                imageSource=
                    "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"
                     +
                    "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"
                     +
                    "7K3xK+Nf/BU7wFqvxY0f4l678Bfhbq91d23w60/w7omk32ra5ZRt5Satqkup2t0oEskbtHbW8cBSNvnllJVh6r/wWCs9f1L/AIJZftA23he1uL3W7nwHq0MEECb5ZFa2dZAq9z5ZfAHJ9D0rzdv2w/Bf7Ef/AARu+GfiDwjcWmtXuoeBtK0L4d6NpgF3c+KdYlsEjsrO2ijBaV2lwX2r8irIzABTXNJ+7WnZtxUFFLdufPt3l7ijFd5PRy5WtutOCaXM5Nt7JQ5Hv0XvXk+0d7OSfoH/AASU/bK179uf9i7SvGPiu0srbxZpuran4b1p7GBobO8urC7ktmuIUYsVSVUV9u47WZlydtfS9fIn7APws0H/AIJBf8EufAXh/wCJmu2un3+jWwufEF0sb3U1/rWoXDTy21vHCrS3c7XExhijiR5ZSqBFJIWu88X/APBSn4beGv2T/Hnxigj8dXvhn4dtcwatayeCdZsdUguII1d4ms7i1juEUB03TPGIowSzuqoxXqxUo03JyafIvecdrqyk0lsnJ6Jd0lujHDxdSygmlNvlT3s23FNt7qK1u3s23o2e/wBFfIf7KX7dHjPXfgnffFr4x3GgeHvAur2dnLouiWHgHxJpviTTr6QSy3Gntb3atPqwjhMHl3NnaxicrcssKpGCfQfB/wDwU7+Bvj25tV0nxwt1a3vhxfFcGonR7+PTJLA/Z+ReNALfzwbu1BtvM+0K1xEpjBcAk4OEnCW8d1vbRvW2mybvtZXTtqKM4ySlHZ7Pvtt82lbdPRq573RXzr4X/wCCr/wE8YzeFk0/xreyDxfrA8PWbyeGtWhjs9TNzJapY37vbKunXLzxSRpDemF5Ch2K1fRVTZ2v0vb5q116q6+9FXV7f11X6P7mFFFFIYUUUUAFFFFABRRRQAUV8L/FH4o/GLQv+C8Hw8+HNr8Vtes/g94s8B3viybwxHomktGbyymjt3gF29m1z5EnmRysBL5gYsFdVKqvQ6Nr3xh+Nv8AwVw8R2nhP4u61p/wF+Fmj2UPirw8NA0qaC/8RTxvKNOgvXtjcqiWrW1xOPMLo08aqyiTEapP2ig19rm+Si3Ft+XNGyte7aXVBV9xzT15eX5uSi0l5pSu+iSbvofY9FFFMAor8w/+CgHjz9qn4GfDTVPFq/GnxZ4H8VfEr4w2fgb4ceD9M0LwzqOm6dpdzfC3hmuXlsJ7iaaS1huLrAuB5fmRowyjgfV/7OPwh+OPww/a58Xp4v8Ait4s+I3whHhbT/7DbxDp2gW922sSXNz9sw2m2drJsihitsCSPaTcthmKnaUf3kFPZba9GoRnZ+dpJf4nbzCt+7k47ta6dVzuF13V03/hVz6NooooAKKKKAPgD/g6O/5QUfHP/uAf+pDpleAf8GVP/KLLx9/2VXUf/TRo9e//APB0d/ygo+Of/cA/9SHTK8A/4Mqf+UWXj7/squo/+mjR6AP1+ooooAKKKKACiiigAooooAKKKKACvwB/4PnP+bXf+5r/APcLX7/V+AP/AAfOf82u/wDc1/8AuFoA/X//AIJO/wDKLL9mn/slXhf/ANNFrVf9pD9hy6+PH7fX7Nnxph8RW+m2vwE/4Sf7RpT2Zlk1n+2NOjs12yhwIvKKbzlW3A4461Y/4JO/8osv2af+yVeF/wD00WtcP+2H+1h43+EX/BVP9jj4XaFqVvbeDfjB/wAJr/wlFo9pFLJe/wBm6RDdWm2RlLxbJXYnYRuBwcigD63r5Y/4Y28Xf8Ptf+Gg/M0f/hA/+FH/APCvPL+0N/aH9pf29/aGfK2bfJ8n+PfndxtxzX1PXzR/w21r/wDw+M/4Zx/snR/+EW/4U1/wsn+08Sf2h9s/tv8As7yPveX5Pl/N93du/ixxQB9L18sfsm/sbeLvgp/wU1/a2+LusSaO3hT42/8ACHf8I8lvcM92n9laVLaXXnoUATMjjZhm3LknHSvY9C/a9+E3ij4vXHw+0z4ofDvUfHto7xz+GrXxJZzavCyAl1a0WQzKVAJIK8Y5ryz9l79trX/jp/wUY/ak+DWo6To9noXwJ/4RP+yL23En2vUP7X0yW8m8/cxT5HQKmxV+UnOTzQB6n+158Dpf2nf2Tvih8NYNRj0ef4h+EtV8Mx38kJmSya9s5bYTFAVLhDJu2gjOMZHWj9kP4HS/sxfsnfC/4az6jHrE/wAPPCWleGZL+OEwpetZWcVsZghLFA5j3bSTjOMnrWP+398WNb+Av7B/xs8deGbmOz8SeC/AWu67pVxJCsyQXdrp088LlGBVgJEU7WBBxgjFH7APxY1v49fsH/BPx14muY7zxJ408BaFruq3EcKwpPd3WnQTzOEUBVBkdjtUADOAMUAfi7/wfOf82u/9zX/7ha/X/wD4JO/8osv2af8AslXhf/00WtfkB/wfOf8ANrv/AHNf/uFr9f8A/gk7/wAosv2af+yVeF//AE0WtAHg/wDwT51a6vP+C5P/AAUMtZrm4ltbP/hW/wBnheQtHBu8Pzltqnhcnk46mvvevgD/AIJ4f8p1/wDgor/3TX/1Hrivv+gD4I/4NhtWutc/4IbfBC6vbm4vLqX+3t808hkkfHiDUgMsck4AA+go/wCC8urXWk/8MZfZbm4tvtP7UvgiCbypCnmxt9u3I2OqnAyDwcVX/wCDXH/lBR8DP+4//wCpDqdH/BfT/myv/s6rwN/7fUAff9fBH/BBrVrrVv8Ahs37Vc3Fz9m/al8bwQ+bIX8qNfsO1Fz0UZOAOBmvvevgD/ggX/zep/2dV45/9saALH/Bzzq11of/AAQ2+N91ZXNxZ3UX9g7JoJDHImfEGmg4YYIyCR9DX3vXwB/wdHf8oKPjn/3AP/Uh0yvv+gD4I/4J86tdXn/Bcn/goZazXNxLa2f/AArf7PC8haODd4fnLbVPC5PJx1NH/BeXVrrSf+GMvstzcW32n9qXwRBN5UhTzY2+3bkbHVTgZB4OKr/8E8P+U6//AAUV/wC6a/8AqPXFH/BfT/myv/s6rwN/7fUAff8AXwR/wQa1a61b/hs37Vc3Fz9m/al8bwQ+bIX8qNfsO1Fz0UZOAOBmvvevgD/ggX/zep/2dV45/wDbGgCx/wAF5dWutJ/4Yy+y3Nxbfaf2pfBEE3lSFPNjb7duRsdVOBkHg4r73r4A/wCC+n/Nlf8A2dV4G/8Ab6vv+gD+YH/gyp/5Sm+Pv+yVaj/6d9Hr+n6v5gf+DKn/AJSm+Pv+yVaj/wCnfR6/p+oAKKKKACiiigDxn9vD9hbwF/wUJ/Z08QfD7x3omj36anYXMGl6ndabDd3Xh27liaNL21MgJjmQkEMhUkAjOCa5v/gnJ8H/AIw/s+/ATwn4A+JFn8JrDTvAWgwaFZy+DpbmX+3GiCot48T21rFY/Ih3W8aThnlLCWMJsf6Kooh7jk4/atf1V0n6pSa+fdRaJ+8o3+ze3zs2vRtJ/Ls2n+dngz9nH49WP/BKf9or4ex/Dq50b4oeNrvxTdWssmvaa7eIrvV9Qu5HmtdkzRQwi0mhEZuZo5TIrq6RhFkk9N+O37NnjD9o74Zfs1+C5/htB4f8EeHvG2n6v4p0e61Kyum0bSdKtLiWwhuERjFJI9zHZJJFbG4jQ7gsjovm19jUUU/c5Wvs+zt/3C1h93bbskFX94pKX2vaf+VVaf3797n5n/tA/sV/Hf4weAv23/D1p4R1DT9T+MusW02j6/D4gsIP7d0K1g0+3i0qzQyytFLLbxakjvdrbxpJcx48xHdovef2K/2ffFHgvxL4m+IXiXwv8S7rxZBoK6FoF18QvEehy65Haphzp8On6EiaNZ2hkigZZkla4mcuJgiRRZ+t6KmMeWHJF2fKo366R5Lp9Jct79HeWnvSvUpczvJXXM5W6atSa/w3SsullroreAf8EvP2ctT/AGXP2GvAnhrxJo9ponja6tH1rxZbwPHLnWb2Rrm9LyRsyysJpGXeGYEIuGKgV7/RRWs5KUm0rLstkuiXktl5GcY2Vm7+fd9W/N7sKKKKgoK+GfFH/BLzVv2eP+CkUH7Qn7Pfg/4M2jeJfDVx4d8W+H9XjOgp573AnTVrS4tLG4b7QxLLPGVjEwVd0m7DL9zUUkrTU1ur/c04tfNNr8VZpNNu8XB7O33ppp/JpP8AB3TafyX+2v8ACX4j3n7Qf7NvxDsvB178XNL+Ft3q1x4h0Lw9Lp2m3EmoXOm/Z7bU7eHU7uKHZCxuV8s3XmILsEGTaxrA/wCCrnjLxD48/YS0T4d6vpFl4e8a/HzxZpHgaLRrXUzfmK1ur9HvN0ixJvKaZDcPMFUomJFDyIA7/adc5f8Awd8I6r8ULDxvdeFfDlz400uyfTrLX5dMhfVLS1clngjuSvmpExJJRWCkk5FNbxT1XMpPzV7tf9vW5bvZbaJRFsrx0ai0n562ff3W722bWu7ZxX7bHh/xhqH7FnxM0b4ZaPHqvja98K32neHdPWeK1SS6kt3ihXfK6RoFLA5Z1A29aj8G+H/DP7AH7DlpZw2tppXhb4R+DjJLHAoSOOGytC8rcA8t5bsTgkliTkmvXKzvF3hDSfiB4W1HQte0vTtb0TV7aSzv9Pv7ZLm1vYHUq8UsTgq6MpIKsCCCQRWVaM3TqKm7Smlr5rm5X8uZl0vZqdPnV4wb0WmkuW6/8lVux+Z37Cv7NPxM/aJ/ZN/Z48H+IfhVq3ws8P6N4ltPi3458Q6pqml3Mni7UhctqdullHaXE85M93JFLNJeLbvHFCIwshYhP1CrK8DeA9D+F/g/TvD3hrRtJ8O6Bo8C2thpmmWkdpZ2MS8LHFFGAiIOyqABWrXVUlHWNJWjdtL1SX4RjGNlZWitL3bwhGdlKo7ysk382398pSeuuu9rWKKKKyNAooooAKKKKACiiigD80v+Csnxl1X9mn/grJ+zJ4n8O6Lea94t8WeD/F3g/wAN2EEcjx32rTnTzaJcFQRHbI7eZLJ/BGjMR8or6L1yTRv+CTf/AATo1rUrvWoLjWdItLnUtT8Salo9/f2+seI712eTUL+PToJbhYJr2UF2SM+XGQowFWvonWfhx4e8ReMtF8R6hoOjX3iHw2twmk6pcWUUt5pa3Cqk4gmZS8QkVVDhCNwUA5wK8E/4KsfCnxp8a/2XbDw74M8LXvjQXHjDw/ea/o1ldWVvdX+j2upwXd3FEbyaC3LMkAXa8qAhmGT0MQjy0o0b7y5W+0ZVG/lbnbk+to3Xuq9tp1HVa2V7d3GCX4qKSW93LX3tPfPh82tP4C0RvEsulT+ImsIDqkmmRSRWT3Xlr5pgSVmkWIvu2hyWC4ySc1zWnftO+BtW/aU1L4QW+ueZ8RdH0CHxPd6T9juB5OnSzNBHP5xTyTmRSuwOXGMlQOa6zwhe6pqXhfT7jW9PtNK1eaBHvLO2vDeQ2spGWRZikZkAPG7YufSuK+Hn7Lvhf4d/Hfxr8TIlvtT8b+O47a0vdTv5VkktLC2UiDT7ZVVVitkdpJNoBZ5JWZ2c7du0nerdqy1vb00S+dn25U1u0YRTVJJO70/NXb+V7dbtdLnyp+3hqNl8df8AgsV+yJ8KDJbTr4KGt/FXVbZpSHQW1sbPT3CgdftM0jDnpE1fQv7dH7TF3+zj8NdDi0LUdB0vxp458Q6d4Y8MzeING1XUNGe9ubmNPLuW0+J3iLRGXyzI0UbSBA0iruI6bVf2PPhHrvxri+JV78LPhzefEaCSOaPxVP4aspNbjeOMRxuLwxmYMsYCKQ+QoAHAxXkf/BRP4afETx/8Uf2fNQ8HeCZvHeg+CvG8viPXbGLU7Kx8qWPTbuDT5pXuZFPkR3VwkkjQrNKgiDJDKwC1nT92FOEv505eac1d+qglHt7qb3dtZq85zj0haPrGLaT8nNt97St0R9QrkKMkE98DFcx8XtJ8aa14NeDwDr/hfw14hMqFL3X9An1yyWMH51NvDeWbliOjecAO6tXS2zSPbRmZEjmKguqOXVWxyASASM98DPoKfVNakrY8A/4Vz+1P/wBFk/Z//wDDN6v/APNPR/wrn9qf/osn7P8A/wCGb1f/AOaevf6KQz8vv+C+f7Pv7TXxI/4JHfGPTNT8efCLxnpi2mn317pWhfDi+0O+lt7bVLO5llW8uNfuYYUhSJpnLwuDHE4+UkOvxx/waP8Axv8AFN7+zF4x+FngP41/BbwX4sn8Y3evL4U8V+B77WtXv4HsLCI3dvLFrFijxZgZTEsbvGYmZm2yIB+8/wAWPhR4b+Onw21rwf4w0Ww8ReGPEVo9jqWm3sfmQXcLjBVh+oIwQQCCCAa/AL/gqV/waSeLPghr03xO/ZE1bVtRh0ucagnhCfUDFrWkujbw+nXhZTNsIyqSMsw2ja8rkCgD9n/+Fc/tT/8ARZP2f/8Awzer/wDzT0f8K5/an/6LJ+z/AP8Ahm9X/wDmnr8T/wDgl/8A8HY/j/8AZk8TRfCz9rrRNd1yz0ib+zpfExsWg8R6I6kIUv7ZgpuAmPmbCzjDFhMxr9+/2f8A9orwL+1V8LdO8a/DnxVovjHwtqq5t9R0y4E0RbAJjcfejkXIDRuFdTwyg8UAeaf8K5/an/6LJ+z/AP8Ahm9X/wDmno/4Vz+1P/0WT9n/AP8ADN6v/wDNPXv9FAHgH/Cuf2p/+iyfs/8A/hm9X/8Amno/4Vz+1P8A9Fk/Z/8A/DN6v/8ANPXv9FAHgH/Cuf2p/wDosn7P/wD4ZvV//mno/wCFc/tT/wDRZP2f/wDwzer/APzT17/RQB4B/wAK5/an/wCiyfs//wDhm9X/APmno/4Vz+1P/wBFk/Z//wDDN6v/APNPXv8ARQB4B/wrn9qf/osn7P8A/wCGb1f/AOaevxA/4PKfDnxT8P8A/DOP/CzPGXw/8W+d/wAJN/Zv/CMeDbzw79kx/ZHm+d9o1S/87dmPbt8rZsfO/cNn9H1fgD/wfOf82u/9zX/7haAP0P8A+CZPgH9pK8/4Jt/s+TaF8WPgfp2iS/DXw4+n2l/8KNUvbq1tzpdsYo5Z08RQpNIqbQ0ixRhiCQiA7R1HxM/YJ+OHxd/aF+GPxR134p/Ae58ZfB/+1f8AhF7tPhNrkUdl/aVstrd7o18UhJd8SKBvB2kZGDXpH/BJ3/lFl+zT/wBkq8L/APpota8H/wCCg2rXVn/wXJ/4J52sNzcRWt5/wsj7RCkhWOfb4fgK7lHDYPIz0NAHvH/Cuf2p/wDosn7P/wD4ZvV//mnr5g/au/4JpftM6v8AFj4kftA+GPjX8Obz4tXXwX1L4ZaXpWj/AA0vtMjmgaeS/ia2mk16Zre/NyVRJ28yJMqTCxHP6PV8Uf8AC6/F3/ERl/wrn/hI9Y/4QP8A4Zw/4ST+wPtLf2f/AGl/wk/2f7Z5WdvneT+734zt46UAfyR/An4VfFHWP2qfDfhXwNpHieD4tw6/DBpdlbwyQapZalHMCpIYBonjkXczPjZsJYgAkf2ifs+/sRaR8CP2rfjR8Y4tZ1PUPFHx1i8PjXbSQRrY2T6RYvZxG2AUOBIrszB2bnGMDivYIPCul2viCfVotNsI9VuY1hmvVt0FxKi9FaTG4qOwJwK+N/2G/jX4u8c/8Fjv26vBuseI9Y1Pwp4F/wCEB/4R7Sbi5aS00X7Xok0115EZOI/NkUO+PvMATQB9Z/GT4T6J8evhD4q8C+JraS88N+NNHu9C1W3jmaF57S6heCZA6kMpMbsNykEZyDmj4N/CfRPgL8IfCvgXwzbSWfhvwXo9poWlW8kzTPBaWsKQQoXYlmIjRRuYknGSc15R/wAFV7uWw/4JeftIzwSyQTwfC3xPJHJGxV42Gk3RDAjkEHnIo/4JUXct/wD8EvP2bp55ZJ55/hb4YkkkkYs8jHSbUliTySTzk0Afj5/wfOf82u/9zX/7ha+Qfj7/AMHJX7SXgj9kT4VfA/wLo1z8DtH8LeANC0g6uIZDr+v28WnQRJewzyoot7edUEsbQpv2sCJmBr6+/wCD5z/m13/ua/8A3C1+kv7LP7Cnwj/by/4I4fs1+GPi14D0Hxppi/CnwyLd7yErd6ezaPagvb3KFZoHPdo3UnocjigC98Gv+Cos3xE+K/jfwR4f/Zp+OGp+Pvh9b6SfGFtDf+DEnsze2pnsmlkfXlEnmQqzLhmKjg7TxXpn/DZHxF/6NO/aA/8ABv4H/wDmhr56/wCCcVmmn/8ABcr/AIKH28eRHAnwzjXJycDw9OBX6CUAfIH7Ln/BUOT9p34E6F45+Fn7Mnxw1zwHrn2j+y72xvvBdtBN5VxLBNtjfX0ZcTRyqcqMlSeQck+Pv/BUOT4D/wDCFf8ACwP2ZPjhpX/CZ+K7Hwt4b+1X3guf7ZrV15n2WCPbr7bJG2SYdtqjBywzXn//AAa4/wDKCj4Gf9x//wBSHU6P+C+n/Nlf/Z1Xgb/2+oA+gP8Ahsj4i/8ARp37QH/g38D/APzQ15/8Av8AgqHJ8eP+E1/4V/8AsyfHDVf+EM8V33hbxJ9lvvBcH2PWrXy/tUEm7X13yLvjy67lORhjivr+vgD/AIIF/wDN6n/Z1Xjn/wBsaAPQP2o/+Cocn7MXwJ13xz8U/wBmT44aH4D0P7P/AGpe3194LuYIfNuIoId0aa+7NmaSJRhTgsDwBkegf8NkfEX/AKNO/aA/8G/gf/5oa+f/APg6O/5QUfHP/uAf+pDplff9AHyB8O/+CocnxB+O3xF8DeHP2ZPjhd+PPh5/Zn/CX2UN94Ljn0/7bbtPZebIdfCyb4VZl2s20DBweKPj7/wVDk+A/wDwhX/CwP2ZPjhpX/CZ+K7Hwt4b+1X3guf7ZrV15n2WCPbr7bJG2SYdtqjBywzXn/8AwTw/5Tr/APBRX/umv/qPXFH/AAX0/wCbK/8As6rwN/7fUAfQH/DZHxF/6NO/aA/8G/gf/wCaGvP/AIBf8FQ5Pjx/wmv/AAr/APZk+OGq/wDCGeK77wt4k+y33guD7HrVr5f2qCTdr675F3x5ddynIwxxX1/XwB/wQL/5vU/7Oq8c/wDtjQB6B8ff+CocnwH/AOEK/wCFgfsyfHDSv+Ez8V2Phbw39qvvBc/2zWrrzPssEe3X22SNskw7bVGDlhmvQP8Ahsj4i/8ARp37QH/g38D/APzQ18//APBfT/myv/s6rwN/7fV9/wBAH8qX/Bop8S9a+Ff/AAUk8bahoXw98YfEq7m+Gt/bvpnhu50qC6gQ6ppTGdm1G9tITGCqqQsjPmRcIVDMv9F3/DZHxF/6NO/aA/8ABv4H/wDmhr8AP+DKn/lKb4+/7JVqP/p30ev6fqAPAP8Ahsj4i/8ARp37QH/g38D/APzQ0f8ADZHxF/6NO/aA/wDBv4H/APmhr3+igDwD/hsj4i/9GnftAf8Ag38D/wDzQ0f8NkfEX/o079oD/wAG/gf/AOaGvf6KAPAP+GyPiL/0ad+0B/4N/A//AM0Nez+APE174y8G6dqmo+HtY8J317EJJtI1WS1kvbBsn93K1rNPAW7/ALuV1561sUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAfIP/AAU+/wCCJHwN/wCCqfhmR/G2hHQ/HMEPlaf4x0ZEg1W2wPkSU423MIxjy5QcAtsMZO6vwN+M37GP7aX/AAbG/Gifx74G1u91T4dzzokviHSYHufD2rxbgEh1SyYnyHOdoL9C5EMxOTX9WFVta0Wz8SaRdafqNpbX9hexNBcW1zEssNxGwwyOjAhlIJBBGCDQB+Y//BJP/g6F+D//AAUB/szwf8QzZfCL4rXOyBLO+uv+JLrkpwo+yXT42OzdIJsN8wVHlOTX6g1+KP8AwVq/4NEPBvxx/tLxv+zTcaf8PvFb7p5/B92xXQdRbgkWzjLWbn5sLhoSSoAhUE18T/sO/wDBej9qD/giH8VI/gz+0T4V8TeKvCGhlIH0LxC5TW9Fg4VXsLxtyzQBR8kbM8TBQI3jGTQB/UNRXhn7CP8AwUf+D3/BSL4XL4q+E/i6z12KJV/tDTJf3GqaPIwz5dzbMd6HOQGGY2KnY7AZr3OgAooooAKKKKACvwB/4PnP+bXf+5r/APcLX7/V+AP/AAfOf82u/wDc1/8AuFoA/X//AIJO/wDKLL9mn/slXhf/ANNFrXf/ABE+PXw++H3x2+HXgbxHq2n2njz4h/2n/wAIhZTWryT6h9it1nvfKkCFY9kLKzbmXcDgZPFcB/wSd/5RZfs0/wDZKvC//pota83/AG1P2bvG/wATf+CuH7E/xD0LQLjUPBvww/4Tr/hKNTSWJY9J+36NDb2m5WYO3mSqyjYrYI5wOaAPsevP/wDhFvhl/wANT/235Hg//hcv/CKfYfO3wf2//YH2zft25877H9q5zjZ5nfdXoFfDH/Cvdf8A+Ilz/hK/7D1j/hFv+GZf7J/tj7FJ/Z/2z/hKvN+zefjy/O8v5/L3btvOMc0Afc9ef/Dvwt8MtJ+O3xF1PwtB4Pj+I+rf2Z/wnEmnPAdWl8u3ZdP+3hD5gxAWEXmAZTO3ivQK+GP2Cvh7r/hz/gtX+33r+o6HrFhoXiT/AIV5/ZGo3FlJFaap5GgzxzeRKwCS+W5CvsJ2sQDg0AfY/wAWfHvh74VfCzxN4o8W3dvYeFPDelXWqazczxNLFb2UELyzyOihiyrGrkgAkgYAPSj4T+PfD3xV+FnhrxR4Su7e/wDCniTSrXVNGuYImiiuLKeFJYJERgpVWjZCAQCAcEDpXnf/AAUU+G+t/GT/AIJ9/HXwh4ZsJNV8SeK/h7r+j6VZRuqPeXdxptxDDEGYhQWkdVyxAGeSBR/wTr+G+t/Bv/gn38CvCHiawk0rxJ4U+HugaPqtlI6u9nd2+m28M0RZSVJWRGXKkg44JFAH4y/8Hzn/ADa7/wBzX/7ha/X/AP4JO/8AKLL9mn/slXhf/wBNFrX5Af8AB85/za7/ANzX/wC4Wv1//wCCTv8Ayiy/Zp/7JV4X/wDTRa0Aeb/sV/s3eN/hl/wVw/bY+Ieu6Bcaf4N+J/8Awgv/AAi+pvLE0erfYNGmt7vaqsXXy5WVTvVck8ZHNfY9eEfAD9uO1+PH7bv7QXwWh8O3Gm3XwE/4Rz7Rqr3glj1n+2LCS8XbEEBi8oJsOWbcTnjpXu9AHxx/wQA/Zu8b/sjf8Ej/AIS/Dz4jaBceF/GXh7+2P7Q0yeWKWS287Wb+4iy0bMh3RSxtwx4bnByKP+Cw/wCzd43/AGjv+GWf+EJ0C417/hA/2hfCfjDXvKlij/s7SbT7X9pum3su5U8xMquWO7hTXqH/AATN/bjtf+CkX7EXgn402Xh248J2vjP7ds0qe8F3Ja/Zb+5szmUIgbcbcv8AdGA+OcZJ+3X+3Ha/sQ/8Kc+1eHbjxF/wt74oaJ8M4fKvBbf2ZJqXn7btso3mLH5JzGNpbd94YoA93r44/wCCPH7N3jf9nH/hqb/hNtAuNB/4Tz9oXxZ4w0HzZYpP7R0m7+yfZrpdjNtV/LfCthht5UV9j14R+wr+3Ha/tvf8Lj+y+Hbjw7/wqH4oa38M5vNvBc/2nJpvkbrtcIvlrJ5wxGdxXb945oA8v/4L/wD7N3jf9rn/AIJH/Fr4efDnQLjxR4y8Q/2P/Z+mQSxRSXPk6zYXEuGkZUG2KKRuWHC8ZOBX2PXhH/BTL9uO1/4Ju/sReNvjTe+HbjxZa+DPsO/SoLwWkl19qv7azGJSjhdpuA/3TkJjjOR7vQB8cfsV/s3eN/hl/wAFcP22PiHrugXGn+Dfif8A8IL/AMIvqbyxNHq32DRpre72qrF18uVlU71XJPGRzR/wWH/Zu8b/ALR3/DLP/CE6Bca9/wAIH+0L4T8Ya95UsUf9naTafa/tN029l3KnmJlVyx3cKa9Q+AH7cdr8eP23f2gvgtD4duNNuvgJ/wAI59o1V7wSx6z/AGxYSXi7YggMXlBNhyzbic8dKP26/wBuO1/Yh/4U59q8O3HiL/hb3xQ0T4Zw+VeC2/syTUvP23bZRvMWPyTmMbS277wxQB7vXxx/wR4/Zu8b/s4/8NTf8JtoFxoP/CeftC+LPGGg+bLFJ/aOk3f2T7NdLsZtqv5b4VsMNvKivsevCP2Ff247X9t7/hcf2Xw7ceHf+FQ/FDW/hnN5t4Ln+05NN8jddrhF8tZPOGIzuK7fvHNAHl//AAWH/Zu8b/tHf8Ms/wDCE6Bca9/wgf7QvhPxhr3lSxR/2dpNp9r+03Tb2XcqeYmVXLHdwpr7Hrwj9uv9uO1/Yh/4U59q8O3HiL/hb3xQ0T4Zw+VeC2/syTUvP23bZRvMWPyTmMbS277wxXu9AH8wP/BlT/ylN8ff9kq1H/076PX9P1fzA/8ABlT/AMpTfH3/AGSrUf8A076PX9P1ABRRRQAUUUUAeUftt2/xT1H9mjxBY/Ba4g0/4laq9pYaTqU620kWiia6hjnv3juAY5RbwNLN5ZVjJ5YUKSwFfHHw0i+PXxS/4KQ/F/4Qab+1H8X5PCfwk8G6Rd3d+fDXg4Xlxr+oebMkW/8AsXy1txbIjeWU37nz5mAM/o8SFBJIAHU18H/8EM9QsvjZa/tHfHS1ktrqL4w/FnVH066ikL+dpWmrHp9nk4HGIZWHtIKiEOerKF9oSl+MYJdlbnc07Xcoq7aSSqcuWkpW3lGP5zfrdU+W3aTa63+nP2INK+Jui/sj/D6D4z6qdZ+Kv9jQy+KLnyLSH/TnBeSPbaAW/wC7LeXmL5W2ZBOcn1SvnrwD+0r4u1//AIKEfFv4e3U/hy++H/w98I6NraCx0O7Gs2t9fNdZtpJRcSJcgRWbSgRW8b/6TGvJXdJy/wALv+Cn3g7Tf2ZdS+L3xF8YeHo/A+sa/q3/AAjFxoPhrXftaaLZ3DQebe2c1v8Aa1mhMbm5lECW8WV+bbiRtpVPayc7JX1stLXlyqNujb+FdUtCFD2a9nvay7393m362W/Z7n1bRXlvwo/bR+G3xv8Ai9rPgTwx4gm1HxLodiuqTQtpV5b213aGUwfaLS6liW3vIlmVo2e2kkVHBViGGKwNC/4KQ/BfxN8VX8H2XjIzX6y3tumpHSL9PD9xNZQtNeQxau0I06WW3jSQyxpcM8flSBgCjgZ3Vk+6b+Sdm/RPRvoxp3vbo0vm1dL1a1XdHuNFfPOkf8FUvgbr9zFb2firWbi9k8Rw+FXs18Jaz9rtb+Y24hE8P2XzILeQ3dqEupVW3f7RHtlO4V7/AKrb/bNLuYvIt7rzYmTyZziKbII2ucN8p6Hg8HoelErqHOlp089E/wAmn6NdwjZy5b/1e35p/cT0V8Af8O5P+rBf2AP/AAof/wAEK9l/Yo/ZT/4UR8RdT1P/AIZn/Zk+C32rTzbf2t8N9U+1ajeZkRvs8q/2Fp+ITt3E+a/zIvyfxLcVd2ZMnZaHof7bdv8AFPUf2aPEFj8FriDT/iVqr2lhpOpTrbSRaKJrqGOe/eO4BjlFvA0s3llWMnlhQpLAV8cfDSL49fFL/gpD8X/hBpv7Ufxfk8J/CTwbpF3d358NeDheXGv6h5syRb/7F8tbcWyI3llN+58+ZgDP6PEhQSSAB1NfB/8AwQz1Cy+Nlr+0d8dLWS2uovjD8WdUfTrqKQv52laasen2eTgcYhlYe0grKEOerKF9oSl+MYJdlbnc07Xcoq7aSS0nLlpKVt5Rj+c363VPlt2k2ut/pz9iDSvibov7I/w+g+M+qnWfir/Y0Mvii58i0h/05wXkj22gFv8Auy3l5i+VtmQTnJ9Urxj4+/8ABQT4U/s0eLpdB8Va5rL6vaWsV9f2uh+GNV8QPo9vK/lwzXv9n204s45W3CNrgxiTY+0tsbHmX7eX/BRa9+AHxz+H/wAH/Bum6v8A8J38QEnvjrNz8Ptf8T6VoWnwKC85ttOjR712laKExx3EYt/PWWZ0XYkutSp7SfPFfE3ZLa+t0vTXTpYzjBQXI38Ku79ujfXXTXrfzPrSivAD+3d4T+Az+DvBfxn8Y+GdM+Kep6av9rDQ9M1FtBXUEs2up4Irl0dIWaKKWSGC4lWeSNQVRs1Hc/8ABU34Fab4O8Oa9qHjS60fTfFGpHSbU6p4e1TT5rScXrWBN7DPbJLYQ/a0aDzrtYovMG3fkgUnbm5U76206u7WnzTXqindK8tOv6n0HRXmH7P/AO2T8Ov2n/EvifRvBmt3t9q3g42x1WzvdGvtKnhjuVdra4RLuGIzW8yxu0c8W+KQKSrsK9Poaa3EmnsFFFFIYUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFeK/tw/wDBPT4Rf8FFvhW/hH4s+D7DxHZRhjY3uPJ1HSJGxmW1uF/eRNwuQDtfaA6sOK9qooA/mD/bs/4N8v2l/wDgjL8UT8Zv2a/FPijxd4V0JmuItV0HMXiLQ4QQzJeWsY23MGAA7Rq0bKrGSKNeD9ff8Ek/+DvXwv8AFv8AsvwN+07BZeC/Ej7Le38a2MRXRtQbhQbyIZNo5PJkXMPJJEKiv2+r8zf+Ctn/AAbH/Bv/AIKJDUvF/gqO0+E3xZud0z6pp1qP7L1uU8/6bargbmOczxbZMsWcS4C0AfpL4c8Sad4x8P2WraRf2Wq6VqUCXNpeWc6z291E4DJJHIpKurAghgSCDxV2v5Sfhb+1D+2t/wAGwPxth8H+LtJur74eXl0zpompSveeF9eTcC82nXa/6iYrydm1lLKZoTgLX71f8Etv+C6PwO/4Kp+H4LTwrrH/AAjHxDig8y/8Ga1KkeoxYALvbt926hBz88fzAYLpHkCgD7MooooAK/AH/g+c/wCbXf8Aua//AHC1+/1fgD/wfOf82u/9zX/7haAP1/8A+CTv/KLL9mn/ALJV4X/9NFrVj4//ALcdr8B/23f2ffgtN4duNSuvj3/wkf2fVUvBFHo39j2Ed426IoTL5ofYMMu0jPPSq/8AwSd/5RZfs0/9kq8L/wDpotar/tIfsOXXx4/b6/Zs+NMPiK3021+An/CT/aNKezMsms/2xp0dmu2UOBF5RTecq24HHHWgD6Hryf8A4bJ8I/8ADcv/AAz55esf8J5/wgv/AAsPzPs6/wBn/wBm/wBof2fjzd+7zvO/g2Y287s8V6xXyx/wxt4u/wCH2v8Aw0H5mj/8IH/wo/8A4V55f2hv7Q/tL+3v7Qz5Wzb5Pk/x787uNuOaAPqevJ/hF+2T4R+Nf7UXxf8AhFo8esL4r+CX9jf8JC9xbqlo/wDato93a+Q4cl8Rod+VXa2AM9a9Yr5Y/ZN/Y28XfBT/AIKa/tbfF3WJNHbwp8bf+EO/4R5Le4Z7tP7K0qW0uvPQoAmZHGzDNuXJOOlAHtf7U3xxi/Zi/Zj+I3xKn06TWIPh54X1PxNJYRzCF71bK0luTCHIYIXEe3cQcZzg9KP2WfjjF+07+zH8OfiVBp0mjwfEPwvpniaOwkmEz2S3tpFciEuAocoJNu4AZxnA6VX/AGvPgdL+07+yd8UPhrBqMejz/EPwlqvhmO/khMyWTXtnLbCYoCpcIZN20EZxjI60fsh/A6X9mL9k74X/AA1n1GPWJ/h54S0rwzJfxwmFL1rKzitjMEJYoHMe7aScZxk9aAPxF/4PnP8Am13/ALmv/wBwtfr/AP8ABJ3/AJRZfs0/9kq8L/8Apota/ID/AIPnP+bXf+5r/wDcLX6//wDBJ3/lFl+zT/2Srwv/AOmi1oA+f/8Agnh/ynX/AOCiv/dNf/UeuK+/6+eP2b/2HLr4D/t9ftJ/GmbxFb6la/Hv/hGPs+lJZmKTRv7H06SzbdKXIl80vvGFXaBjnrX0PQB8Af8ABrj/AMoKPgZ/3H//AFIdTo/4L6f82V/9nVeBv/b6veP+CUn7Dl1/wTd/YF8BfBa98RW/iy68Gf2hv1WCzNpHdfatRurwYiLuV2i4CfeOSmeM4B/wUK/Ycuv23v8AhR32XxFb+Hf+FQ/FrQfiZN5tmbn+049N+0brRcOvltJ5wxIdwXb905oA+h6+AP8AggX/AM3qf9nVeOf/AGxr7/r54/4J6/sOXX7EP/C8ftXiK38Rf8Le+LWvfEyHyrM239mR6l9n22jZdvMaPyTmQbQ277oxQB4P/wAHR3/KCj45/wDcA/8AUh0yvv8Ar54/4Kt/sOXX/BSL9gXx78FrLxFb+E7rxn/Z+zVZ7M3cdr9l1G1vDmIOhbcLcp94YL55xg/Q9AHwB/wTw/5Tr/8ABRX/ALpr/wCo9cUf8F9P+bK/+zqvA3/t9XvH7N/7Dl18B/2+v2k/jTN4it9Stfj3/wAIx9n0pLMxSaN/Y+nSWbbpS5Evml94wq7QMc9aP+ChX7Dl1+29/wAKO+y+Irfw7/wqH4taD8TJvNszc/2nHpv2jdaLh18tpPOGJDuC7funNAH0PXwB/wAEC/8Am9T/ALOq8c/+2Nff9fPH/BPX9hy6/Yh/4Xj9q8RW/iL/AIW98Wte+JkPlWZtv7Mj1L7PttGy7eY0fknMg2ht33RigDwf/gvp/wA2V/8AZ1Xgb/2+r7/r54/4KFfsOXX7b3/CjvsviK38O/8ACofi1oPxMm82zNz/AGnHpv2jdaLh18tpPOGJDuC7funNfQ9AH8wP/BlT/wApTfH3/ZKtR/8ATvo9f0/V/MD/AMGVP/KU3x9/2SrUf/Tvo9f0/UAFFFFABRRRQBQ8UeF9M8ceGtQ0XWtOsNY0fV7aSyvrC9t0uLa9gkUpJFLG4KujKSrKwIIJBGDXLfBD9mT4bfsy6VfWPw2+Hvgf4fWOqSrPeW/hrQrXSYruRRtV5Ft0QOwHALAkCu4ooWl2uoPVJPofIH7N/wAMvi38Nfjx+1ZrF34Ca2vfHmuzaz4b8QXGs2RtNWhg0qzs9LtoYUkkmQqYJjM1ysIQsmwSh28ry3Sv2RPi34t/4IneEvgJc/Debw7rt6vh/wAL+IrS71vTprtrBry1bXdRlMMr2w3j7bIiRzTvIjoWUSO0K/ojRShFRioNXSVNa9VTulf/ABJ2n3XZ6jk25OS0d5vTpz9v8L1h281ofC3xN/Zg+L15+2r8bb7wh4WTRdH8S/CG08CfD/xZHqtra6d4WMcN/JIgt1ke5843k1kVC26xCKBm80MixyZ3/BL79ibxX8Mpfho3j3wl8V7aP4SeGE0fQ18e+KfDzQeHbprdYJv7KsdAjMN1HIjTI93qcoukVYwit507D76oqqTcG2tb9/Wcr+t6ktevXd3molNJbW7eShH7rQjp5Hzn/wAE5v2fNe+D3hv4m+JfGnh208P+Ofif8Qdb8RaiiTQ3Er2f2p7fTQ8sTupAsYbdgu47PMYEK24V9GUUUlpGMFtFKK9IpRX4JDespS/mcpP1k3J/iwooooAoeKPC+meOPDWoaLrWnWGsaPq9tJZX1he26XFtewSKUkiljcFXRlJVlYEEEgjBrlvgh+zJ8Nv2ZdKvrH4bfD3wP8PrHVJVnvLfw1oVrpMV3Io2q8i26IHYDgFgSBXcUULS7XUHqkn0Pjr9jnwF8V/2dPjz8ZdK1H4V3mst8Rvibe+Kv+E+m17TLfSLjRpUto7eEokj6h9qtraPyEha0ETND/x8Ir7hZ+F+oD46f8FnPiTrcHl3WjfBL4f6d4OjmBLpDq2qXLaheRqcbVkW2g0/eAc4kTPUY+vCMiuc+F3we8I/A/wy2i+CvC3hzwfoz3Ml41hommw6fbNPId0kpjiVV3ueWbGSepNFL3HT/wCnceVPr8Hs1fp8Dl03aYVPeVT+/Lmf/gftHb1kl5WueLfH/wDZ81741/8ABQX4G69e+HbS/wDh38MNM1zXpNQuJoWEOvTLa2tiqwl/MLLBJeuHEZVSB8wYqD5n/wAFf/Fl5rvjf9nP4a6b4D1T4ny+JvHyeJtU8MafdWUE+o6Xotu91If9NmgtnCXUli+yWZA20AZJFfbFcz4t+Cng3x/448O+J9e8JeGNb8SeEHlk0LVr/S4Lm+0VpQFka1mdS8BcABjGV3ADOaUVZ07fZkpX6tqXMn2dmorpeKSvfUbs1O/2ouPycXFr8W+tm27dDxf9jL4F+L7b48fFv41+PtKPhXxB8VH03TtN8MNdwXc+gaPpkcyWq3UsBaFruWS5uZZBDJLHGJI0WWTYWP0dRRVN6Jdlb+v18ybauT3f/DL7lZLyQUUUUhhRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFAHI/HP4B+Cv2mvhlqXgz4geGNG8X+FtXTZdabqdss8MnowB5V16q6kMpwQQRmvwH/4Klf8ABpT4x+AOvzfE/wDZG1fWdXtdKn/tFPCct8Ytd0d0O8Pp13lTPsI+VGKzDauHmY1/RNRQB/Oh/wAEu/8Ag7L8d/s2eI4fhb+1xo2ua5ZaRP8A2bJ4o+xND4i0R0IQpqFswU3ATB3OAs4wxImY1+/nwB/aI8DftT/C3TfGvw68U6N4w8Lasm621HTbgTRE4BKMPvRyLkBo3AdTwwB4r5u/4Kh/8EPfgd/wVT8OSz+MNF/4R7x7DB5Wn+MtGjSLU4MA7Em423MIP/LOXJAzsaMndX4HfGD9kT9tL/g2I+Nk3jfwXrV5qHw9ubhI31/S4Xu/DetR7sJDqdm2fIlIO0b8EF2EMxOWoA/q1r8Af+D5z/m13/ua/wD3C19g/wDBJP8A4OgPg7/wUF/szwf4/ey+EfxXudkCWN/dAaPrcpwB9junwFdmPEE2HywVGlOTXx9/wfOf82u/9zX/AO4WgD9f/wDgk7/yiy/Zp/7JV4X/APTRa1w/7Yf7WHjf4Rf8FU/2OPhdoWpW9t4N+MH/AAmv/CUWj2kUsl7/AGbpEN1abZGUvFsldidhG4HByK7j/gk7/wAosv2af+yVeF//AE0Wtdv8TP2T/BHxd/aF+GPxR13Tbi58ZfB/+1f+EXu0u5Yo7L+0rZbW73RqwSXfEigbwdpGRg0AekV80f8ADbWv/wDD4z/hnH+ydH/4Rb/hTX/Cyf7TxJ/aH2z+2/7O8j73l+T5fzfd3bv4scV9L14v/wAMS6B/w8P/AOGjv7W1j/hKf+Fdf8K2/szMf9n/AGP+0/7R8/7vmed5ny/e27f4c80Ae0V80fsvftta/wDHT/gox+1J8GtR0nR7PQvgT/wif9kXtuJPteof2vpkt5N5+5inyOgVNir8pOcnmvpevF/gj+xLoHwL/bA+OPxl07VtYvNd+O39g/2vZXBj+yaf/ZFk9nD5G1Q/zo5Z97N8wGMDigC5+398WNb+Av7B/wAbPHXhm5js/EngvwFruu6VcSQrMkF3a6dPPC5RgVYCRFO1gQcYIxR+wD8WNb+PX7B/wT8deJrmO88SeNPAWha7qtxHCsKT3d1p0E8zhFAVQZHY7VAAzgDFdx8ZPhPonx6+EPirwL4mtpLzw3400e70LVbeOZoXntLqF4JkDqQykxuw3KQRnIOaPg38J9E+Avwh8K+BfDNtJZ+G/Bej2mhaVbyTNM8FpawpBChdiWYiNFG5iScZJzQB+Ev/AAfOf82u/wDc1/8AuFr9f/8Agk7/AMosv2af+yVeF/8A00WtfkB/wfOf82u/9zX/AO4Wv1//AOCTv/KLL9mn/slXhf8A9NFrQBw/7Hn7WHjf4u/8FU/2x/hdrupW9z4N+D//AAhX/CL2iWkUUll/aWkTXV3ukVQ8u+VFI3k7QMDAr63r4A/4J4f8p1/+Civ/AHTX/wBR64r7/oA+SP8AghT+1h43/bi/4JWfC34o/EbUrfV/GXij+1v7Qu4LSK0jl+z6ve2sWI41VFxFDGOAMkZPJNH/AAVo/aw8b/spf8Mz/wDCE6lb6d/wsn49eFvAmvebaRXH2rSb77V9phXep8tm8pMOuGXHBGa83/4Ncf8AlBR8DP8AuP8A/qQ6nR/wX0/5sr/7Oq8Df+31AH3/AF8kf8El/wBrDxv+1b/w0x/wm2pW+o/8K2+PXinwJoPlWkVv9l0mx+y/ZoW2KPMZfNfLtlmzyTivrevgD/ggX/zep/2dV45/9saAPSP+C637WHjf9h3/AIJWfFP4o/DnUrfSPGXhf+yf7Pu57SK7ji+0avZWsuY5FZGzFNIOQcE5HIFfW9fAH/B0d/ygo+Of/cA/9SHTK+/6APkj9jz9rDxv8Xf+Cqf7Y/wu13Ure58G/B//AIQr/hF7RLSKKSy/tLSJrq73SKoeXfKikbydoGBgUf8ABWj9rDxv+yl/wzP/AMITqVvp3/Cyfj14W8Ca95tpFcfatJvvtX2mFd6ny2bykw64ZccEZrzf/gnh/wAp1/8Agor/AN01/wDUeuKP+C+n/Nlf/Z1Xgb/2+oA+/wCvkj/gkv8AtYeN/wBq3/hpj/hNtSt9R/4Vt8evFPgTQfKtIrf7LpNj9l+zQtsUeYy+a+XbLNnknFfW9fAH/BAv/m9T/s6rxz/7Y0Aekf8ABWj9rDxv+yl/wzP/AMITqVvp3/Cyfj14W8Ca95tpFcfatJvvtX2mFd6ny2bykw64ZccEZr63r4A/4L6f82V/9nVeBv8A2+r7/oA/mB/4Mqf+Upvj7/slWo/+nfR6/p+r+YH/AIMqf+Upvj7/ALJVqP8A6d9Hr+n6gAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKqa9oNj4q0S70zVLK01LTdQhe3urS6hWaC5icFWjdGBVlYEggggg1booA/En/AIK1/wDBob4S+M39p+OP2Zbiw8C+J33XE/g29kK6HqLdSLWTlrNzzhDuhJKgeSoJr8XP26bL9q3S/hj4d8B/H7SfiInh/wCCOo3Gj6TJ4jsncaPLfxwsbZLwg+dE8dijRASOgVD5ZCtz/azRQB/Kl+z1/wAHdf7SX7NfwC8D/DnQvBPwPu9E8AeH7Dw3p89/o+qSXU1vZ20dvE8rJqKI0hSNSxVFBJOFA4rl/i7/AMHS/wC0P8Zv2ovhB8WNR8NfCey134L/ANs/2RYWOnanHpupf2paJazfbI2v2eXy0QNHsdNrkk7hxX9atFAH8wP/ABGrftT/APQg/s//APgj1f8A+WdcR/xF2ftV/wDDSX/Cxfsnw38n/hGv+Ec/4RP7Jqn/AAje77V9o/tH7N9v3/bsfufN8zb5Xy7M/NX9WtFAH8wP/Eat+1P/ANCD+z//AOCPV/8A5Z1xHgT/AIO7P2q/A/xn8e+MWtPhvrcfjr+z9mgapaapNovhr7JA0J/s6EX6tB5+7zJtzvvkUEbRxX9WtFAH8q3x6/4O9v2l/wBoX4GeNPAGreDfgnp2l+ONCvvD95d6bpOqw3trDd27wSSQO2osqyqshKsysAwBII4o+Av/AAd7ftL/ALPXwM8F+ANJ8G/BPUdL8D6FY+H7O71LSdVmvbqG0t0gjknddRVWlZYwWZVUFiSABxX9VNFAH8YX/BVv/gtX8U/+Cwf/AAgX/CzNA+H+h/8ACu/7Q/s3/hGLG8tvP+2/ZfN877Rcz7sfZI9u3bjc+d2Rj6Q/Z6/4O6/2kv2a/gF4H+HOheCfgfd6J4A8P2HhvT57/R9Ukuprezto7eJ5WTUURpCkaliqKCScKBxX9VtFAH8nXwz/AODqb44fCL9oX4nfFHQvhf8AAe28ZfGD+yv+Eou30/XJY73+zbZrW02xtqhSLZE7A7ANxOTk16P/AMRq37U//Qg/s/8A/gj1f/5Z1/T9RQB/J1+yP/wdTfHD9h39nrw/8Lvhz8L/AID6R4N8L/af7PtJ9P1y7ki+0XMt1LmSTVGdsyzSHknAOBwBR+0b/wAHU3xw/at/4QP/AITb4X/AfUf+FbeMNP8AHeg+Vp+uW/2XVrHzPs0zbNUHmKvmvlGyrZ5BxX9YtFAH8wP/ABGrftT/APQg/s//APgj1f8A+Wdecfs5f8HU3xw/ZS/4Tz/hCfhf8B9O/wCFk+MNQ8d695un65cfatWvvL+0zLv1Q+WreUmEXCrjgDNf1i0UAfydftcf8HU3xw/bi/Z68QfC74jfC/4D6v4N8UfZv7QtINP1y0kl+z3MV1FiSPVFdcSwxngjIGDwTXo//Eat+1P/ANCD+z//AOCPV/8A5Z1/T9RQB/J18M/+Dqb44fCL9oX4nfFHQvhf8B7bxl8YP7K/4Si7fT9cljvf7NtmtbTbG2qFItkTsDsA3E5OTR+0b/wdTfHD9q3/AIQP/hNvhf8AAfUf+FbeMNP8d6D5Wn65b/ZdWsfM+zTNs1QeYq+a+UbKtnkHFf1i0UAfzA/8Rq37U/8A0IP7P/8A4I9X/wDlnXnH7OX/AAdTfHD9lL/hPP8AhCfhf8B9O/4WT4w1Dx3r3m6frlx9q1a+8v7TMu/VD5at5SYRcKuOAM1/WLRQB/J1+0b/AMHU3xw/at/4QP8A4Tb4X/AfUf8AhW3jDT/Heg+Vp+uW/wBl1ax8z7NM2zVB5ir5r5Rsq2eQcV6P/wARq37U/wD0IP7P/wD4I9X/APlnX9P1FAH8yH/BlF4W1K5/4KTfEjW47C7fR7L4aXdjPeiImCGeXVNLeKJn6B3WGZgOpEbHsa/pvoooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAK8Q0n9sr+1P+Ci2sfAH/hG/L/snwBb+Of7e/tDPm+bfvZ/Zfs/l8Y2b/M805zjYOte318QeE/8AlYt8Y/8AZAdO/wDT9PRT1r04PZ89/lTnJfik/l2FV0w9Sa3XJb51YRf4Nr59z6j8SfEDxxpn7Qvhnw5p3w+/tTwFqmmXd1q/jD+3beD+xLqMr5Fr9hYedP5wLHzEIVNnPUV3lfGvx8/5Tmfs6/8AZNPGf/pTo9fCP7KX/BHj4R/tw/suftI/ED4oz+OvE2t6P8R/Hq+E4f8AhJ7y2sfBMkN9cO0+n20brCs0soSSQypIrmKP5Rht2LrKNFVZbKFSb72hV5LLo3rptorO7vJ6qnes6a3coRXa8qfP6paa76u60tFftxRX4d/tc/H74r/tD/8ABNr/AIJ/eErrw5rXxg0z412IHjXQv+E6i8IzePZ7awja20641WY8ee7PK6bt87QbVIcqw9T/AOCSnwG+JH7Kv7Qnxc8Iv8CJf2aPhRrPw7fV18BXPxgs/HW3VUmaIalbxCZrm2SaIyRyMVMbtbRjeCqqNcT+5dZP/l3zpdLuEXJ3va10rK3M72uknczoP2saMl/y8UH3spy5Ft2ervyq2zb0P1xor8yv+Ddb/gnN8MfCX7HHwU/aCmtfEeufF3WfBwsP7d1TX7ydbHT2YounwWokW2S2RYxtXyiwJLFiTmv01rfEUfZVHTb1Tafyb/r1utVq8cPW9rTVRLRpNfNL+vTXR6IooorE2CiiigAooooAqeIPEFh4T0K91TVb200zTNNge6u7y7mWGC1hRSzySOxCoiqCSxIAAJNfIPgH/g4K/Y5+KHx3tPhtofxz8PXniu/1FtJtVbT7+HT7q5VioSO/kgWzcOy4R1mKyFlCFiy557/g4LmbXf2OPBHgm8vJ7Hwr8UPil4W8IeKZYXeNzpNzfqbhN6kFQ3lqpPcMR3r6/HwJ8ED4Z6Z4LPg/wu3g/RBarp+htpcDadYi1dJLbyoCvlp5Txxsm1RsZFIwQDRR969SXwqXLbrpGMm7+k1bTV31XV1fdShH4nHmv0V3KK066xd9VZWte7t80/tO/wDBfD9kz9jb43638OPiR8V/+Ec8Z+HDEuo6d/wjGs3n2cyxJMn723tJImzHIh+VzjODggitr4/f8Frf2ZP2XvhZ8OfGvjr4l/2H4Z+LWmtq/hS8/wCEd1W5/tW1VYnMnlw2zyRfLNEdsqo3zdODj498G/tJfH/4D/8ABWz9s6D4L/s0/wDC+7PUda8MSapdf8LD03wt/Y8i6JGI49l2jGbeCxymAuzB6iul/wCClXxr+Lvgf/gol+xP4v8AB3wR/wCE5+KNx4V8WS3PgD/hMbLTPskstjY/aYv7SlUwP9n3P8wXEmz5cZFRSblTpN7zSb0fWnKdlHfdJc2yWr0aKqJKpKK2V+q/mUb823W/Lu3otUz7g/ZA/wCCg3wa/b4+GV/4v+EXj3SvGeh6VI0V+0EU1vdaew3YE9tOiTxbgjFd8Y3hSV3DmvEPgx/wcRfsafH74s6R4H8M/G7SpfEmvXP2Oxh1HRNV0mCabBIj+0XdrFArMRtUM4LuVVcsyg+Sf8Efri//AGj/ANoT9rD9oDxRpWifDnx74ivbfwTr/wANLAyT3XhWbSrd187ULpooku7i48zcksCGLyVQCR23rH8V/s5/E/xn/wAFAP8AgiR8NP2SfAX7OPxhvfEWsSQwJ8R/EXhuOx8D6HDHqstxNqltqTyEyyRR70VUVJGYyBNzAI+1Nc1RRtuqbsmvt815c3w8qSUtrWer0uZuyUuZ2SlJc1n9m2jj8V7txsnduOi95I/oCoqto2nnSdItbUzS3BtoUiMshy8u1QNzH1OMmrNTJJNpCg24pyVn2CiiikUFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABXy9+2h/wRh/Zq/4KF/FK08a/GD4bf8ACX+JrHTo9JgvP+Eg1XT9lskkkix+Xa3MUZw0sh3Fd3zdcAY+oaKlxTabW235FKTSaT33/M8M/Z6/4Jr/AAU/ZUbwGfAPgv8AsE/DLTdT0jw1/wATe/uv7NtdSuUub2P99O/m+ZMitul3suMKVGRXX/Cj9lLwB8EPh14l8J+F9B/svw/4v1PUdY1e1+3XM/2u61B2ku5N8kjOnmM7HajKq5+UKK9Eoqqnv359bpp36pvmafk3q+71epMfdd46O6fzSsn6paJ9tDw3xv8A8E1vgb8Sv2QtH+A3iH4eaTrXwp8PWsNppei3s9xO2mrEpWN4blpDcxyqrMolWUSYdhuwxBzP2N/+CU/7Pv8AwT+0XxHY/CL4bab4RTxank6tcfbbu/vLuLbt8o3NzLLMsXfy1cIGJbG4k19C0US95ylLVy+Lz9e/zCOiiltHby9O3yOP+APwE8J/su/Brw98PvAuk/2H4Q8K2gsdLsPtU1z9lhBJC+ZM7yNyTy7E+9dhRRVSk5Nyk7tkxiopRirJBRRRUlBRRRQAUUUUAeZfti/sj+C/26f2b/FHws+IFlPe+GPFVsIJzbSCK5s5FYPFcQOQQk0Uio6Eqy7lAZWUlT8p/DL/AIJY/tN+GfEHhzSPEv7eXxE8SfC3w7fW8w0O38E6dpuv39rbSCSC2n11ZHuZMlI1nkZSZ08xWCiQ4++aKIe5Lnj5PyutnbZtd7XCfvR5JefrrvZ7q/qeIfs6/sa/8KC/ap+O/wATf+Ek/tb/AIXZqGk3/wDZv9n+R/Y32GxFps83zG87zMb87I9vTDdaPjD+xr/wtj9uH4N/Gb/hJPsH/CpNP12w/sf+z/N/tX+04IYt/n+YvleV5WceW+/djK4yfb6KIuzi19lWXkuXk/8ASXb8d9RNXTT6u79b8356nz14e/YKj8D/ALfHxA+NegeJk0yy+KfhK00DxP4c/ssOt/f2jOLXU1uBKu10gkaFkMbbhg71xg73/BPT9kT/AIYN/Y08C/CP/hIf+Er/AOEKtJbX+1fsH2H7Z5lxLNu8nzJNmPN243t93PfA9nooh7seSO2n4OTX3Ocvk7bJJOXvT53v/wACK/KK+6+7dyiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKAP/Z",
                fileName=
                    "U:/Documents/PHD/Articles/1st paper/gfx/tank_wall.jpg")}));
      end HeatTransferTank_no_PCM;

      model HeatTransferTank_PCM

        import SI = Modelica.SIunits;

      /******************** Connectors *****************************/
        Ports.HeatFlow2 heatFlow
          annotation (Placement(transformation(extent={{-26,80},{-6,100}})));

      /****************** General parameters *******************/
      parameter Integer tank=4 "Tank type" annotation (Dialog(group= " Heat Transfer"),choices(
      choice=1 "Tyep 1 Steel",
      choice=2 "Tyep 2 Aluminium",
      choice=3 "Tyep 3 Aluminium liner",
      choice=4 "Type 4 Plastic liner"));

        inner parameter Boolean Charging = true
          "If true, tank is charging with given heat transfer coefficient" annotation(Dialog(group= " Heat Transfer"));
        parameter SI.CoefficientOfHeatTransfer  h_charging=150 "Charging, Heat transfer 
  coefficient inside HSS tank 150 w/m2K - 500w/m2K  
    according to monde (0 is adiabatic fuelling)"                                                             annotation(Dialog(group= " Heat Transfer"));

        parameter SI.CoefficientOfHeatTransfer  h_discharging=-1 "Discharging heat transfer 
  coefficient, if <0 then Daney relation
     is used, if >0 then the given number is used "                                                           annotation(Dialog(group= " Heat Transfer"));
        inner parameter SI.CoefficientOfHeatTransfer  h_o=8
          "Heat transfer coefficient outside the tanks (typical natural convection)"
                                                                                     annotation(Dialog(group= " Heat Transfer"));
       // inner parameter Real  T_amb=273;

       inner parameter SI.Length xLiner=0.003 "Thickness of liner" annotation(Dialog(group= "Geometry"));
       inner parameter SI.Length xCFRP=0.022 "Thickness of wrapping/tank" annotation(Dialog(group= "Geometry"));

      outer SI.Length d;

      protected
       inner parameter SI.Length x1=xLiner/(t1-0.5);
       inner parameter SI.Length x2=xCFRP/(t2-0.5);
       inner constant Real t1=5.5;
       inner constant Real t2=10.5;
       inner Real y1;
       inner Real y2;

      SI.Temperature Delta_T1;
      SI.Temperature Delta_T11;
      Real Delta_T22;
       Real Delta_T33;

      public
       inner SI.CoefficientOfHeatTransfer  h_i = h_charging;
       //(if Charging == true then h_charging else h_discharging);

      inner parameter Boolean Adiabatic_Wall = true;

      parameter Integer N = 100 "number of PCM volumes";
      inner parameter SI.Length t_PCM = 1e-2 "total PCM-layer thickness";
      final inner parameter SI.Length dx_PCM =  t_PCM/N; // discretization step
      inner SI.Length   Lc = (2*(log(((d/2+t_PCM))/(d/2)))^(4/3))/((d/2)^(-3/5) + (d/2+t_PCM)^(-3/5))^(5/3);
      final inner parameter SI.Temperature T_amb = 293.15;

      inner parameter Boolean Convection = true;

      inner parameter Real Go = 1
          "multiplying coefficient for PCM thermal conductivity";
      inner parameter Real G1 = 1 "multiplying coefficient for PCM density";
      inner parameter Real G2 = 1
          "multiplying coefficient for PCM specific heat capacity";
      inner parameter Real G3 = 1
          "multiplying coefficient for PCM latent heat of phase change";

      /************************************************************************************/

      parameter Boolean Manual_Area_inner = false;
      parameter Boolean Manual_Area_PCM = false;
      inner Real  G_A_inner = (if Manual_Area_inner == false then 1 else G_Inner)
          "multiplying coefficient for heat transfer area at the H2/PCM interface";
      inner Real  G_A_PCM = (if Manual_Area_PCM == false then 1 else G_PCM)
          "multiplying coefficient for heat transfer area within the PCM layer";

      parameter Real G_Inner = 2
          "multiplying coefficient for heat transfer area at the H2/PCM interface";
      parameter Real G_PCM = 2
          "multiplying coefficient for heat transfer area within the PCM layer";

      /************************************************************************************/

      /****************** equations *******************/

        PCM_NEW_with_Nat_Conv_Coupling_Area.PCM_inner pCM_inner(each Delta_T1=Delta_T1, each Delta_T22=Delta_T22)
          annotation (Placement(transformation(extent={{6,50},{26,70}})));
        PCM_NEW_with_Nat_Conv_Coupling_Area.PCM pCM[N](each Delta_T1=Delta_T1, each Delta_T22=Delta_T22)
          annotation (Placement(transformation(extent={{-64,30},{-44,50}})));

        WallPieces_PCM_Area_middle.Liner5Pieces wall_liner
          annotation (Placement(transformation(extent={{-4,14},{16,34}})));

        WallPieces_PCM_Area_middle.Tank10Pieces wall_CFRP
          annotation (Placement(transformation(extent={{-40,-16},{-20,4}})));

        WallPieces_PCM_Area_middle.OuterWallCell outer_wall
          annotation (Placement(transformation(extent={{-10,-38},{10,-18}})));

      equation
      //Deciding which calls to make in lookup tables in wallpieces for tank properties
       if tank==1 then
         y1=1;
         y2=1;
       elseif tank==2 then
         y1=2;
         y2=2;
       elseif tank==3 then
         y1=2;
         y2=4;
       elseif tank==4 then
         y1=3;
         y2=4;
       else
         y1=3;
         y2=4;
         end if;

      Delta_T1 = pCM_inner.T_PCM - pCM[N].B.T;
      Delta_T11 = pCM_inner.A.T - pCM[N].B.T;

      Delta_T22 = pCM_inner.T_PCM - 329.15;

      Delta_T33 = pCM_inner.T_PCM - 293.15;

        connect(heatFlow, pCM_inner.A) annotation (Line(
            points={{-16,90},{-16,60.1},{6.5,60.1}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(pCM_inner.B, pCM[1].A);
        for i in 1:N-1 loop
         connect(pCM[i].B, pCM[i+1].A);

        end for;
        connect(pCM[N].B, wall_liner.heatFlow);
          connect(wall_liner.heatFlow1, wall_CFRP.heatFlow);
        connect(wall_CFRP.heatFlow1, outer_wall.portA);

        annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
                  -100},{100,100}}),
                            graphics), Icon(coordinateSystem(preserveAspectRatio=
                  false, extent={{-100,-100},{100,100}}),
                                            graphics={Bitmap(
                extent={{-90,-80},{96,74}},
                imageSource=
                    "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"
                     +
                    "jNKlP4OVx06J78vbS69JSWl7rWFZwk6q1ldPXq1qr/NRs9bWWjtZ91+xx8DYP2Zf2Tfhr8PLdESPwV4Z0/Rm2psDvBbpG749WdWY8nknk9a9Joorqr1pVakqst5Nv79TloUlSpxpLaKS+4KKKKyNT4A/4Ojv+UFHxz/7gH/qQ6ZXgH/BlT/yiy8ff9lV1H/00aPXv/8AwdHf8oKPjn/3AP8A1IdMrwD/AIMqf+UWXj7/ALKrqP8A6aNHoA/X6iiigAooooAKKKKACiiigAooooAK/AH/AIPnP+bXf+5r/wDcLX7/AFfgD/wfOf8ANrv/AHNf/uFoA/X/AP4JO/8AKLL9mn/slXhf/wBNFrXz/wD8FD/+U6//AATq/wC6lf8AqPW9fQH/AASd/wCUWX7NP/ZKvC//AKaLWu/+In7UfgT4UfHb4dfDPX9d+weNviz/AGn/AMIppv2K4l/tX+zrdbm8/epG0UXlwsrfvXTdnC7jxQB6BXxB/wALY8Vf8RI//CC/8JN4g/4Qn/hmr+3f+Ee/tGb+yv7R/wCEo8j7b9m3eV9o8n935u3fs+XOOK+368//AOLWf8NT/wDNP/8Ahdv/AAin/Tn/AMJV/wAI99s/8Cv7P+1/9sfO/wBugD0CviD9hL4seKvF/wDwWc/by8Lat4m8Qap4Y8H/APCvv7B0i71Gaew0T7Toc8tz9lgZjHB5sgDyeWF3sAWyea+368/+Hf8Awqz/AIXt8Rf+ET/4V/8A8LN/4ln/AAnn9kfY/wC3f+Pdv7O/tTyv3/8Ax77/ACPtH/LPds+XNAHAf8FYv+UWX7S3/ZKvFH/pouqP+CTv/KLL9mn/ALJV4X/9NFrXr/xZ+KWhfA74WeJvGvim+/svwx4P0q61vV73yZJ/slnbQvNPL5catI+2NGbaisxxgAnAo+E/xS0L44/Czw1418LX39qeGPGGlWut6Re+TJB9rs7mFJoJfLkVZE3RurbXVWGcEA5FAH4Q/wDB85/za7/3Nf8A7ha/X/8A4JO/8osv2af+yVeF/wD00WtfkB/wfOf82u/9zX/7ha/X/wD4JO/8osv2af8AslXhf/00WtAHn/7Gn7Lnjv4Uf8FYv20PiZr+hfYPBPxZ/wCEH/4RTUvttvL/AGr/AGdo81tefukkaWLy5mVf3qJuzldw5r6/rwD9nj9uf/hff7dn7RfwU/4Rb+yf+FA/8I1/xOf7S8/+3f7Y0+S9/wBR5S+R5OzZ/rJN+c/JjFe/0AfIH/BBX9lzx3+xd/wSe+FPwz+Jehf8I1428Nf2v/aWm/bbe8+zefrF9cxfvbeSSJt0M0bfK5xuwcEEA/4K7/sueO/2oP8AhmD/AIQXQv7c/wCFd/tAeFfG3iH/AE23tv7P0ey+1/abr99InmbPNT93Hukbd8qNg49A/wCCXH7c/wDw8o/YT8DfGv8A4Rb/AIQv/hNPt/8AxJv7S/tH7H9l1C5sv9f5UW/d9n3/AOrXG/HOMk/b6/bn/wCGHP8AhSv/ABS3/CUf8Lg+Kuh/DL/kJfYv7I/tLz/9O/1UnneV5P8Aqvk37v8AWLjkA9/r5A/4JEfsueO/2X/+Gn/+E60L+w/+FiftAeK/G3h7/Tbe5/tDR737J9muv3Mj+Xv8p/3cm2RdvzIuRn6/rwD9gX9uf/huP/hdX/FLf8Iv/wAKf+KuufDL/kJfbf7X/s3yP9O/1Ufk+b53+q+fZt/1jZ4APP8A/gvV+y547/bR/wCCT3xW+Gfw00L/AISXxt4l/sj+zdN+229n9p8jWLG5l/e3EkcS7YYZG+ZxnbgZJAP1/XgH/BUf9uf/AIdr/sJ+OfjX/wAIt/wmn/CF/YP+JN/aX9nfbPtWoW1l/r/Kl2bftG//AFbZ2Y4zke/0AfIH7Gn7Lnjv4Uf8FYv20PiZr+hfYPBPxZ/4Qf8A4RTUvttvL/av9naPNbXn7pJGli8uZlX96ibs5XcOaP8Agrv+y547/ag/4Zg/4QXQv7c/4V3+0B4V8beIf9Nt7b+z9Hsvtf2m6/fSJ5mzzU/dx7pG3fKjYOPQP2eP25/+F9/t2ftF/BT/AIRb+yf+FA/8I1/xOf7S8/8At3+2NPkvf9R5S+R5OzZ/rJN+c/JjFH7fX7c//DDn/Clf+KW/4Sj/AIXB8VdD+GX/ACEvsX9kf2l5/wDp3+qk87yvJ/1Xyb93+sXHIB7/AF8gf8EiP2XPHf7L/wDw0/8A8J1oX9h/8LE/aA8V+NvD3+m29z/aGj3v2T7NdfuZH8vf5T/u5Nsi7fmRcjP1/XgH7Av7c/8Aw3H/AMLq/wCKW/4Rf/hT/wAVdc+GX/IS+2/2v/Zvkf6d/qo/J83zv9V8+zb/AKxs8AHn/wDwV3/Zc8d/tQf8Mwf8ILoX9uf8K7/aA8K+NvEP+m29t/Z+j2X2v7TdfvpE8zZ5qfu490jbvlRsHH1/XgH7fX7c/wDww5/wpX/ilv8AhKP+FwfFXQ/hl/yEvsX9kf2l5/8Ap3+qk87yvJ/1Xyb93+sXHPv9AH8wP/BlT/ylN8ff9kq1H/076PX9P1fzA/8ABlT/AMpTfH3/AGSrUf8A076PX9P1ABRRRQAVQ8VaPceIfDGo2Fpql/od1fW0lvDqNisLXVg7KVWaITRyRGRCQyiSN0yBuVhkG/RSlFSTi+o07O6PBP8Agnv+wDof/BN/4K3vgTwv4y8e+MdCn1K41aH/AISqexuLm1nndpZ9sttawM4kkZnPm7yCeCBxX5S+MvEXwV/aK/Z6+GHij4pN4L8QfHD4gfGexk+I3iPWDb3158LbS21a4uRoks8g3abGbexFnDZ5jMzSSuEkLOx/dWinFuNWFVbwcLekJRdl5NRSfktLa3maUqc6b+3zX9ZKSb83eTfr8rVNB1qHxJodlqNul3Hb38CXMSXVrLazqrqGAkhlVZI3weUdVZTkEAgirdFFN2voCvbUKKKKQwooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigD4A/wCDo7/lBR8c/wDuAf8AqQ6ZXgH/AAZU/wDKLLx9/wBlV1H/ANNGj17/AP8AB0d/ygo+Of8A3AP/AFIdMrwD/gyp/wCUWXj7/squo/8Apo0egD9fqKKKACiiigDwH9r39l7xD+0j8WfhzPcePfE3gz4WeD49T1PxPZ+GvFWpeHdS8QXLRRR2cMlxZSQuLSMG5lk/fK29IQAVLkfCn/BPT9j7xx+3J/wT3sfiZpXxh+Omi+JfGfxGl1nQbu9+MPieSHTvCkGuqhsliN3JHMz2EE6q08bOzyjdIo5X7E/4LR/tc6T+xX/wTO+LXi/UNVsdN1O60G60XQY55gkl9qd1C8MEUS5DO4LGQhckJFI3RCRgfBD47+Bf+CcH/BGX4ceMbGC/+IHgH4feCtKW7vPBxs7szwiKOO4vU3zxRuiyF3fY5f72FYgijCu06lRfYdP0vKTk7+aUIryhJp6O5VeLkqdO1+fn9bRSVl3u6jffmjG2yPsSivGfEP7c3g3Rf2ldW+Elvb6zq3jfSPAb/EKW0tEgWN7AXBt0hEksqKtw7g7VcqmOWdRXU/svftAWP7VX7PXhD4j6XoniLw7pfjTTY9VsrDXYIoNQhgk5jaVIpJUUum1xtdvldehyA4puPMtv+DKP5wkvkzPnV0u/+UZflKL+Z3tFfJ/xZ/4K6+EfhpqTz6b8Ovit448Iw+MLbwDJ4t8P2mmHSDrc12ln9kiW5voLq4Edw4jeW3t5IVcOgkLRuq/WFKPvQ51t/wABP8mmu6aa0ZUtJOD3X+bX5pp9mmnqgr8Af+D5z/m13/ua/wD3C1+/1fgD/wAHzn/Nrv8A3Nf/ALhaAP1//wCCTv8Ayiy/Zp/7JV4X/wDTRa15/wDtl/sueO/iv/wVi/Yv+JmgaF9v8E/Cb/hOP+Er1L7bbxf2V/aOjw21n+6eRZZfMmVl/dI+3GW2jmvQP+CTv/KLL9mn/slXhf8A9NFrVj4//tx2vwH/AG3f2ffgtN4duNSuvj3/AMJH9n1VLwRR6N/Y9hHeNuiKEy+aH2DDLtIzz0oA93r4Y/4V7r//ABEuf8JX/Yesf8It/wAMy/2T/bH2KT+z/tn/AAlXm/ZvPx5fneX8/l7t23nGOa+568n/AOGyfCP/AA3L/wAM+eXrH/Cef8IL/wALD8z7Ov8AZ/8AZv8AaH9n483fu87zv4NmNvO7PFAHrFfDH7BXw91/w5/wWr/b71/UdD1iw0LxJ/wrz+yNRuLKSK01TyNBnjm8iVgEl8tyFfYTtYgHBr7nryf4Rftk+EfjX+1F8X/hFo8esL4r+CX9jf8ACQvcW6paP/ato93a+Q4cl8Rod+VXa2AM9aAK/wDwUJ+Fuu/HH9gX44eCvC1j/anifxh8P9e0TSLLzo4Ptd5c6dcQwReZIyxpukdV3OyqM5JAyaP+Ce3wt134HfsC/A/wV4psf7L8T+D/AIf6Domr2XnRz/ZLy2063hni8yNmjfbIjLuRmU4yCRg10n7U3xxi/Zi/Zj+I3xKn06TWIPh54X1PxNJYRzCF71bK0luTCHIYIXEe3cQcZzg9KP2WfjjF+07+zH8OfiVBp0mjwfEPwvpniaOwkmEz2S3tpFciEuAocoJNu4AZxnA6UAfh7/wfOf8ANrv/AHNf/uFr9f8A/gk7/wAosv2af+yVeF//AE0WtfkB/wAHzn/Nrv8A3Nf/ALha/X//AIJO/wDKLL9mn/slXhf/ANNFrQB8/wD/AATw/wCU6/8AwUV/7pr/AOo9cV9/188fs3/sOXXwH/b6/aT+NM3iK31K1+Pf/CMfZ9KSzMUmjf2Pp0lm26UuRL5pfeMKu0DHPWvoegD4A/4Ncf8AlBR8DP8AuP8A/qQ6nR/wX0/5sr/7Oq8Df+31e8f8EpP2HLr/AIJu/sC+Avgte+IrfxZdeDP7Q36rBZm0juvtWo3V4MRF3K7RcBPvHJTPGcA/4KFfsOXX7b3/AAo77L4it/Dv/Cofi1oPxMm82zNz/acem/aN1ouHXy2k84YkO4Lt+6c0AfQ9fAH/AAQL/wCb1P8As6rxz/7Y19/188f8E9f2HLr9iH/heP2rxFb+Iv8Ahb3xa174mQ+VZm2/syPUvs+20bLt5jR+ScyDaG3fdGKAPB/+Do7/AJQUfHP/ALgH/qQ6ZX3/AF88f8FW/wBhy6/4KRfsC+PfgtZeIrfwndeM/wCz9mqz2Zu47X7LqNreHMQdC24W5T7wwXzzjB+h6APgD/gnh/ynX/4KK/8AdNf/AFHrij/gvp/zZX/2dV4G/wDb6veP2b/2HLr4D/t9ftJ/GmbxFb6la/Hv/hGPs+lJZmKTRv7H06SzbdKXIl80vvGFXaBjnrR/wUK/Ycuv23v+FHfZfEVv4d/4VD8WtB+Jk3m2Zuf7Tj037RutFw6+W0nnDEh3Bdv3TmgD6Hr4A/4IF/8AN6n/AGdV45/9sa+/6+eP+Cev7Dl1+xD/AMLx+1eIrfxF/wALe+LWvfEyHyrM239mR6l9n22jZdvMaPyTmQbQ277oxQB4P/wX0/5sr/7Oq8Df+31ff9fPH/BQr9hy6/be/wCFHfZfEVv4d/4VD8WtB+Jk3m2Zuf7Tj037RutFw6+W0nnDEh3Bdv3TmvoegD+YH/gyp/5Sm+Pv+yVaj/6d9Hr+n6v5gf8Agyp/5Sm+Pv8AslWo/wDp30ev6fqACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooA+AP8Ag6O/5QUfHP8A7gH/AKkOmV4B/wAGVP8Ayiy8ff8AZVdR/wDTRo9e/wD/AAdHf8oKPjn/ANwD/wBSHTK8A/4Mqf8AlFl4+/7KrqP/AKaNHoA/X6iiigAooooAK8w/bY+Atz+1L+x78UPhvZXFpaX3jrwtqOh2k91nyIJ7i2kjjeTCsdiuyk4UnAOBmvT6Kzq01UhKnLZpr7zSjVlTqRqR3TT+4+D4v+CUfjvV/Gfwm1nWPHujPfjQdX0f4t39rBNHdeI4b7+ymNlp7DBhtwmmRWfmSsZRbDcP9Ibzl9g/bB/aq1/9kTWbfUbXwlcXPwe8F+DdU13xje6fpLteWRjNvDptrp7GSO2LsftTyiT93DFAryPAhUv9IUVrXk6qcXs+bT/Fzad7Jybir2T23d8qMY07aXtb8Lfmlb7uyR+Pf7Duk6z8WvEH7MXwA0Hxh8I/iV4H+CPiCfxr4n1r4d6+fENrcRwC9l06fVJ0gWCxvJb2WMrZxzTtKYriUyFI8N+wlFFVKd1bu3J+cnZN/ckrbaaWWiXK3PnbvokvRNu33yb+et3dsr8Af+D5z/m13/ua/wD3C1+/1fgD/wAHzn/Nrv8A3Nf/ALhago/X/wD4JO/8osv2af8AslXhf/00WtV/2kP2HLr48ft9fs2fGmHxFb6ba/AT/hJ/tGlPZmWTWf7Y06OzXbKHAi8opvOVbcDjjrVj/gk7/wAosv2af+yVeF//AE0WtcP+2H+1h43+EX/BVP8AY4+F2halb23g34wf8Jr/AMJRaPaRSyXv9m6RDdWm2RlLxbJXYnYRuBwcigD63r5Y/wCGNvF3/D7X/hoPzNH/AOED/wCFH/8ACvPL+0N/aH9pf29/aGfK2bfJ8n+PfndxtxzX1PXzR/w21r//AA+M/wCGcf7J0f8A4Rb/AIU1/wALJ/tPEn9ofbP7b/s7yPveX5Pl/N93du/ixxQB9L18sfsm/sbeLvgp/wAFNf2tvi7rEmjt4U+Nv/CHf8I8lvcM92n9laVLaXXnoUATMjjZhm3LknHSvqevmj9l79trX/jp/wAFGP2pPg1qOk6PZ6F8Cf8AhE/7IvbcSfa9Q/tfTJbybz9zFPkdAqbFX5Sc5PNAHqf7XnwOl/ad/ZO+KHw1g1GPR5/iH4S1XwzHfyQmZLJr2zlthMUBUuEMm7aCM4xkdaP2Q/gdL+zF+yd8L/hrPqMesT/DzwlpXhmS/jhMKXrWVnFbGYISxQOY920k4zjJ61j/ALf3xY1v4C/sH/Gzx14ZuY7PxJ4L8Ba7rulXEkKzJBd2unTzwuUYFWAkRTtYEHGCMUfsA/FjW/j1+wf8E/HXia5jvPEnjTwFoWu6rcRwrCk93dadBPM4RQFUGR2O1QAM4AxQB+Lv/B85/wA2u/8Ac1/+4Wv1/wD+CTv/ACiy/Zp/7JV4X/8ATRa1+QH/AAfOf82u/wDc1/8AuFr9f/8Agk7/AMosv2af+yVeF/8A00WtAHD/ALHn7WHjf4u/8FU/2x/hdrupW9z4N+D/APwhX/CL2iWkUUll/aWkTXV3ukVQ8u+VFI3k7QMDAr63r4A/4J4f8p1/+Civ/dNf/UeuK+/6APkj/ghT+1h43/bi/wCCVnwt+KPxG1K31fxl4o/tb+0LuC0itI5fs+r3trFiONVRcRQxjgDJGTyTR/wVo/aw8b/spf8ADM//AAhOpW+nf8LJ+PXhbwJr3m2kVx9q0m++1faYV3qfLZvKTDrhlxwRmvN/+DXH/lBR8DP+4/8A+pDqdH/BfT/myv8A7Oq8Df8At9QB9/18kf8ABJf9rDxv+1b/AMNMf8JtqVvqP/Ctvj14p8CaD5VpFb/ZdJsfsv2aFtijzGXzXy7ZZs8k4r63r4A/4IF/83qf9nVeOf8A2xoA9I/4LrftYeN/2Hf+CVnxT+KPw51K30jxl4X/ALJ/s+7ntIruOL7Rq9lay5jkVkbMU0g5BwTkcgV9b18Af8HR3/KCj45/9wD/ANSHTK+/6APkj9jz9rDxv8Xf+Cqf7Y/wu13Ure58G/B//hCv+EXtEtIopLL+0tImurvdIqh5d8qKRvJ2gYGBR/wVo/aw8b/spf8ADM//AAhOpW+nf8LJ+PXhbwJr3m2kVx9q0m++1faYV3qfLZvKTDrhlxwRmvN/+CeH/Kdf/gor/wB01/8AUeuKP+C+n/Nlf/Z1Xgb/ANvqAPv+vkj/AIJL/tYeN/2rf+GmP+E21K31H/hW3x68U+BNB8q0it/suk2P2X7NC2xR5jL5r5dss2eScV9b18Af8EC/+b1P+zqvHP8A7Y0Aekf8FaP2sPG/7KX/AAzP/wAITqVvp3/Cyfj14W8Ca95tpFcfatJvvtX2mFd6ny2bykw64ZccEZr63r4A/wCC+n/Nlf8A2dV4G/8Ab6vv+gD+YH/gyp/5Sm+Pv+yVaj/6d9Hr+n6v5gf+DKn/AJSm+Pv+yVaj/wCnfR6/eD44/tQ/EfwP/wAFPfgp8JfD0vgu98E/EDQdb1rxFBcaLdvrGkQ6ekYSeO7W7WAJNPcwxhXtyVMb/M28bCPvTjTW7v8AgnJ/gmD0i5Pp/mkvvbS9dNz6eooooAKK+f8A4oftU+JfF37TE3wa+Eljod34p0CytdX8Y+Idbjln0nwfaTuTbwm3ieOS7vblI5SkKzQrGg86STHlxTfP37RX/BY74lfB7wB468d6D8CvCviH4f8Ahj4gj4c6Jeah8RZ9L1XxhqH21LB5LS0GkzRiJbszRkvcAkW0rAEAZUGpNJfa28/ejH8ZSSXd3tohyTim3st/L3XL/wBJi2+y3P0Bor5y/Zx/a2+J/wAQf2ufF/wn+IXws8J+EJfCnhbT/Ex1nw943n8QWlx9tubmCG2KzabZOj/6JcMThhhVAzuyPo2qadlLv+ja/NfrsTf3nHqrfik1+DTCiiikMKKK8V/bb/am139mLwt4Ni8IeCYviH41+IHim08LaJoc2rnSYZJJY5p5riW5EE5jhgt7eeZz5THbGQOSAU3ql3aXzbsvxfy6jSbTfZNv0Su/uSPaqK/PPVP+Cv8A8bdP0v8AaLvofgB8L9Qsf2ZRt8TXFr8Xbxo9SmWxF7LBYltAHmSRxkK4m8kBztBPJH3p8Pdev/FXgHQ9U1XTBouqalp8F1eacLj7QLCZ41Z4fM2rv2MSu7auducDOKqKvHnW1ov5SV0/mlf7r7omT5Zcr3vJfOLSkvVNq/8AwGbFFfP3xL/ar8Sfsy/tG+HtE+I1hocvw0+JerpofhfxVpaSwSaJqkqgW+l6nBI8gb7QyyCG8iZUMhWF4Ijsll+gaUdYqa229Gt0/NXXyaavFptvSXK/X5dGvJ2fzTT1TSKKKKACiiigAooooAKKKKACiiigAooooA8Y/wCCiPx31/8AZd/Yg+J/xI8MX3hyw1zwL4futctW13TJ9RsbhrdDJ9neKCeCTMu3y1ZZBsaRWKuFKHsP2afFXifx3+zt4F1zxtbadZeMNa0Cxv8AWrawt5Le2tryWBJJo445HkdVV2ZQGdjxyTXyl/wcDajqXif9iLQfhboTQf2/8cvH/h7wPZLOrPbsJb1LmbzlR0dofJtZA4V1JViuRnNZ3wcsPFXgL/gtbB4F0r4j/EHxVo2j/CRtd8ex6zrU11pl3qVzqCQ2D29ju+y6fII4JzstY4VaPGQ5LsTDPnlKP80nFeThTdSXyaaXqlZO7sYj3Yxkuiu/NTqRpx+akpfJvXY+9a8x/a3/AGptF/ZD+ELeJ9WstQ1m+vr620XQtE05Va+8Q6pdSCK1sYAxCh5JCMsxCood2IVSR6dXxL/wUC1CbXv+Cq/7D3hW6fzdCm1bxZ4iltSAVkvbHRwtrIc/88/tUxHuQe1JRcpxpp2u/wAFrK3nyp2umr2urFJ2jKVr2UmuzaTaT62bte1nbZo67x1+1l8cfhT4o+GPgR/h38PvHXxV+KM2q6pJptl4hu9E0TwZpFnHCzNc37Wt3LdsktxbwecltAJnmUrBGARXi/gn/gsp8Z/HfwJ1L4k6f+zz8Pr3wrpnxET4bqbP4sXMt3qd22qw6X9rtEbQ1jktvtEwwzyRuVRjsHGftH9rn41W/wCzf+yx8R/H91MsEPgzw1qGs72XcA0Fu8i8YOcsoGMd68d/4IofA64/Z8/4JX/BPQb+ORNXvPDkWu6p5sflyteagzX03mD++HuCpJ67aug06kpSXuw5LrvzSdknvbkpyTbbfNLmvsiKqtCEU9Zc2un2Vq+1+acGkly2i1ZdfqWiio7ppUtZDAkckwQmNHcorNjgFgCQM98HHoahuyuUld2JKK8A/wCFjftT/wDRG/2f/wDw8mr/APzMV2XwU8V/GXXvEVzF8RvAXwy8K6SluWt7nw349vvEFzLNuXCPDPo9iqJt3HeJGOQBs5yKSuS3Y5L9rD9q/wAZfCD4y/DT4dfDj4e6V8Q/GHxEGpXkiar4kk0HT9C06xjiM15cTx2d2+0y3FvCqrESXmXHAOPmXwT/AMFlPjP47+BOpfEnT/2efh9e+FdM+IifDdTZ/Fi5lu9Tu21WHS/tdojaGsclt9omGGeSNyqMdg4z9o/tc/Gq3/Zv/ZY+I/j+6mWCHwZ4a1DWd7LuAaC3eReMHOWUDGO9eO/8EUPgdcfs+f8ABK/4J6DfxyJq954ci13VPNj8uVrzUGa+m8wf3w9wVJPXbSw+s5Slqoct135pNpLsuWnJPd3lzbJRTr/DCMdHLm1/wrV69eacGtLWi0/P6lrzn9pbXviZ4P8AA6a58MtI8N+KtR0Znur7w3qbyWs/iGAIf9Hs7wP5drc5wyNNFLHIVEbGEOZ4/RqKUk2tHYaaT1PzF/4L6/tG+F/2tf8Ag22+KXxB8HXFxPoPiKDQpIkuYTBdWkqeJNOimtp4zny5oZUkikTJ2vGwyRyfO/8Agyp/5RZePv8Asquo/wDpo0evIP2/rl9F/wCCN3/BTDwrbxmHRPDvxzgl06IKoSD7XrGhXUyLjnb50kj4I48zgnt6/wD8GVP/ACiy8ff9lV1H/wBNGj1V1KEKqVlOMJW7c8VK3yva/UmzjOdNu/LKUb9+WTjf52ufr9RRRSGFFFFABRRRQAUUUUAFFFFABX4A/wDB85/za7/3Nf8A7ha/f6vwB/4PnP8Am13/ALmv/wBwtAH6/wD/AASd/wCUWX7NP/ZKvC//AKaLWu3+Jn7J/gj4u/tC/DH4o67ptxc+Mvg//av/AAi92l3LFHZf2lbLa3e6NWCS74kUDeDtIyMGuI/4JO/8osv2af8AslXhf/00WteD/wDBQbVrqz/4Lk/8E87WG5uIrW8/4WR9ohSQrHPt8PwFdyjhsHkZ6GgD73rxf/hiXQP+Hh//AA0d/a2sf8JT/wAK6/4Vt/ZmY/7P+x/2n/aPn/d8zzvM+X723b/DnmvaK+KP+F1+Lv8AiIy/4Vz/AMJHrH/CB/8ADOH/AAkn9gfaW/s/+0v+En+z/bPKzt87yf3e/GdvHSgD7Xrxf4I/sS6B8C/2wPjj8ZdO1bWLzXfjt/YP9r2VwY/smn/2RZPZw+RtUP8AOjln3s3zAYwOK9or4o/Yb+Nfi7xz/wAFjv26vBuseI9Y1Pwp4F/4QH/hHtJuLlpLTRfteiTTXXkRk4j82RQ74+8wBNAH1n8ZPhPonx6+EPirwL4mtpLzw3400e70LVbeOZoXntLqF4JkDqQykxuw3KQRnIOaPg38J9E+Avwh8K+BfDNtJZ+G/Bej2mhaVbyTNM8FpawpBChdiWYiNFG5iScZJzXlH/BVe7lsP+CXn7SM8EskE8Hwt8TyRyRsVeNhpN0QwI5BB5yKP+CVF3Lf/wDBLz9m6eeWSeef4W+GJJJJGLPIx0m1JYk8kk85NAH4+f8AB85/za7/ANzX/wC4Wv1//wCCTv8Ayiy/Zp/7JV4X/wDTRa1+QH/B85/za7/3Nf8A7ha/X/8A4JO/8osv2af+yVeF/wD00WtAHb/DP9k/wR8Iv2hfid8UdC024tvGXxg/sr/hKLt7uWWO9/s22a1tNsbMUi2ROwOwDcTk5NekV8Ef8E+dWurz/guT/wAFDLWa5uJbWz/4Vv8AZ4XkLRwbvD85bap4XJ5OOpr73oA83/ZH/ZP8EfsO/s9eH/hd8OdNuNI8G+F/tP8AZ9pPdy3ckX2i5lupcySMztmWaQ8k4BwOAKP2jf2T/BH7Vv8Awgf/AAm2m3Go/wDCtvGGn+O9B8q7lt/surWPmfZpm2MPMVfNfKNlWzyDivlD/g2G1a61z/ght8ELq9ubi8upf7e3zTyGSR8eINSAyxyTgAD6Cj/gvLq11pP/AAxl9lubi2+0/tS+CIJvKkKebG327cjY6qcDIPBxQB9715v+zl+yf4I/ZS/4Tz/hCdNuNO/4WT4w1Dx3r3m3ctx9q1a+8v7TMu9j5at5SYRcKuOAM16RXwR/wQa1a61b/hs37Vc3Fz9m/al8bwQ+bIX8qNfsO1Fz0UZOAOBmgD6v/a4/ZP8ABH7cX7PXiD4XfEbTbjV/Bvij7N/aFpBdy2kkv2e5iuosSRsrriWGM8EZAweCa9Ir4I/4OedWutD/AOCG3xvurK5uLO6i/sHZNBIY5Ez4g00HDDBGQSPoa+96APN/hn+yf4I+EX7QvxO+KOhabcW3jL4wf2V/wlF293LLHe/2bbNa2m2NmKRbInYHYBuJycmj9o39k/wR+1b/AMIH/wAJtptxqP8Awrbxhp/jvQfKu5bf7Lq1j5n2aZtjDzFXzXyjZVs8g4r5Q/4J86tdXn/Bcn/goZazXNxLa2f/AArf7PC8haODd4fnLbVPC5PJx1NH/BeXVrrSf+GMvstzcW32n9qXwRBN5UhTzY2+3bkbHVTgZB4OKAPvevN/2cv2T/BH7KX/AAnn/CE6bcad/wALJ8Yah4717zbuW4+1atfeX9pmXex8tW8pMIuFXHAGa9Ir4I/4INatdat/w2b9qubi5+zftS+N4IfNkL+VGv2Hai56KMnAHAzQB9X/ALRv7J/gj9q3/hA/+E20241H/hW3jDT/AB3oPlXctv8AZdWsfM+zTNsYeYq+a+UbKtnkHFekV8Ef8F5dWutJ/wCGMvstzcW32n9qXwRBN5UhTzY2+3bkbHVTgZB4OK+96AP5gf8Agyp/5Sm+Pv8AslWo/wDp30ev1g1HwZe/tX/8Fv8A42apP4r1nwt4A+CPw20jwjrlxpV5JYX91JfPNqk8EF7FIslpGYxA0ssGyfMMQSWMBt35P/8ABlT/AMpTfH3/AGSrUf8A076PX9Cfw0/YM8F/DOx+NcEV54i1Rfj3rV5rXiV727QSRtc2kdm1vbPFHG0cKQxgJks6lid54xjUU7ucFrGMuW+zk1ypPy5ZTv6K91oaQ5XHkk7KTinbdRT5m153jFfNtao8v/4IX+NfF/xN/wCCangnxP4w13xF4hl8R3mqajo1zr15Je6mukPqFx9gS4uJCZJ3FuI/3jsxYFfmIxX13X50ftjfAK0/4J1fsK/Dr4dWXxT8f23w51zxd4b8B+J/FniDxDBp3/CKeFEMgaJJbSG2trRZAi2r3QjSd/tYMs7FIymb/wAEwrL4D6B/wVU+PE3wr8MeEfBVqfD3h7w34d0rw7ocdmuq2MdvPqFxrWy3QBbO4a4toku5Qsc7QRbZHLxg9d4VKslTb5YtxTesvdhGWvdvmjd3bb55aqN3ytypwXOkpNKVltaU3Gy7JNSsnZJckb3kkvT/APgjfqE3jLxz+114o1F/tOtah8eNd0mW4IG42unw2lpaRf7scSAAerMe9UP+CqEH/C7v24v2M/g8rGW2vvHl18QtWgEPmIbXQ7J5YzJxgI1xPCvX7xXritvwLY/8O3v21/ivqHiOOWy+C3x51a38V2niMQ/8S7wr4g8iO0vLbUZFXFtFdCK3kiuJSIjL5kbOrtEJPRfFX/BOfSvFn/BQbRP2jZfiR8TLfxR4f0c+HrLQoZdLOgxac5DzW/lPYtP+9lHmNJ5/m5wFdUVUGVF/7tLpBU0/KVOCS/8AKkYvzh7y0cb7VdXiUt5udv8ADUk//bJO3RSVnsz6Jor4E/bflvvgP/wVx+EnirwT4cN544+L3w68ReBbe4gsi8b3cN3plza3F464/c2sbXMzmQj91HIqZdkRvPv+CY/7MFp8ZbP4e+DfEOj6pqumfsc/EDxWsHiTXbU/bNZ1X+07+GziilZR5gS2lS6uHixEZjaBclJUiKP7yMZLre67KM3GT+S5ZJO13NRTvqFZOmnLdaJebceZab2clKLavZRlJq2h+nF00qWshgSOSYITGjuUVmxwCwBIGe+Dj0NeCf8ACxv2p/8Aojf7P/8A4eTV/wD5mK9/opW1uO+h5l8FPFfxl17xFcxfEbwF8MvCukpblre58N+Pb7xBcyzblwjwz6PYqibdx3iRjkAbOcj0fUdQh0nT57q4kEVvbRtLK56IqjJP4AVNXE/tH/BRP2j/AIFeKfAc3iTxJ4StfFunyaXdapoD20epW8Eo2yrC9xDNGjPGWQsYyyhyUKuFZZrufs37JLmtp69L/wDA6eZVGMede0fu318l1sfKX/BADSJvEv7FviL4r3qudR+PPxA8Q+PXkkh8uRree+eG1Huv2e3jZf8AZcV9x185+C/gXd/8E0/+Cbet+E/h3qXir4gz/DDwlqEnhWHXRaSX0zW9rI9rZ5tbaBHAdVUExlzu+ZmNfnH8K5P2cr344fsS+ILbUPCnjXxhNNe+NvG3xOjgj1nV/FHiCPSQF0f7ZGrXF1e/a7+OVbCPc8CW0KiFAIwN1yTrewpfBBQiu/LZxjp1ajB31V5cqV+bTCcpQputU+KbqSt0uvelruk3NW00jzSlbld/t7/g4GsI5v8AgkH8Zr/G298PWFnrenzKAXtru1v7aeCVc9CsiKfzHevrbwNrE3iHwVo+oXKeXcX1jDcSrgDa7xqxHBI6k9Ca+Qf+Ci8Cf8FE57f9mXwWZtW0q812xuPirrlqpOneGdJtLiG8k0x7jBQ6ldlYY1tkLSRxSNLKqIYzJ9oQxLbxLGgCogCqB2ArOl/DlLpJ6L0Wsl/iuo/9w+1jSr8UYreKd/m1ZPzjZu3RTVt2OooooAKKKKACiiigAooooAKKKKACiiigDyz44fsi+G/j/wDGz4T+Otcvtcj1D4O6peazo1lazRLZXdzcWj2pe5Vo2d/LSRmTY6YY5O4cVy95+wH4ds/2wfEXxs0/xT8QtO1rxTZaZb65oOn6nDb6XrLaZ5rWbOwhF2m0yndFFcxwTAASxSKWDe91wv7Tf7PWh/tX/AHxZ8OPEs+qW2heMNPk068l06cQ3MaN3QsrIeQMpIjxuMq6OjMpiaajeno9WvVq1/W1rX7LsitJPlns0k/RPmt6X1Px78JL8K/2+PBn7PDeJLHw548/aR+K3xhg1Txrqlxbx3Wr+BbbT7q71CXQ2lYtNp0UdrYiCOzXYHAllKZd5D+kH/BSr9nzxP44X4W/FfwBpba549+A3igeJbTSImjSfX9NmgktdT0+F3womltZWaMMyq0sUasQDuHUfCP9hKx8B/GbTPH/AIq+IfxC+K3ifw5p0ml+HZvFH9lW9t4ahm4uDa2umWNlbiSZQiNLJHJIEjCIyKXVvdq1dowUaelpc67J+6kkv5bQV073bkrtauG+apKc1o48r80+dt37++0ttk0lsvmr9on4P+AP+Cx/7GGreCofHfjfw94P8RXKWniBPD4t9N1qF4Sskml3sV7ayy2kiuYzLC0ccwKqrEIzK/Pf8FHPFXij9iH/AIJO+Km8KeK/EV3rPhrSrHRF8V3sNsdR0u0mu4LSfU3FtFDCHtraWSXckSqvkhiMBjX1tXkP7Un7H9h+09rfgbWv+Ev8Y+BvEnw61KfVND1bw/8A2fNLbyzW0lrJug1C1u7VyYpXAcw+YmTsdQzhs5pWaitJOLku9vu0s3pfZvW+pdOTTi5PWN7fPXz1dlrbotLKx8X/ALC3wV+BHif/AIK6tr/wO8P+FG8N/Cr4Ux2114p0eOKX/hLNR1a+IW9kvkLPqMgh0+4V7qV3Z5JZfmc7iP0ury39mb9k3Rf2Zv8AhJr+HW/E3jHxd43vxqPiPxR4juIZtU1qVF8uFXEEUNvFFDFiOOG3hijRQTt3M7N6lWrfuRp78t/vlKUn0WicmlotEtDGMfflPvb/AMljGK76tRu/NhRRX5kfte/8HW37PH7Fn7S/jL4V+KfBvxov/EPgjUG02+uNK0nTJbKWQKrExNJfxuVww5ZFPtUGh9o/t7/sVaP/AMFCP2adZ+FXiPxZ408JeGvEbxDVZfDE9pBeX8CNv+zNJc284WJnCM2xVc7Au7Yzq3kf/BUT9l201L/gjF8VfAt7d6l42ufCPgW6vdM1LV7e0fUJLrT4GuLafFvDDEsytEoDRRoeOBknPyH/AMRq37LH/Qg/tAf+CPSP/lnR/wARq37LH/Qg/tAf+CPSP/lnWVSEvZzjTfK5W131SfK/le6NqVVRqwnNXUenk2rr52Sfkl2Oo1v4s+MPiV+2R8PPiSnhjWdZh/aT+EuteBvCPhfUdPkOnw28Vxpc0N/foVPkwTLc393KZwCbSO3jCCfMT/Z2i3nw6/4JNfsaeB/Alkmo6jB4a0xND8M6BpdqbrXfF97HEXaK0tUy0s8rCSVyMRxhpJJGjiRnX4J/4jVv2WP+hB/aA/8ABHpH/wAs6P8AiNW/ZY/6EH9oD/wR6R/8s66JyTi4QXKm+m9uacktdLp1J62s29U0kjlpQ5bOWrS+V+WMe97WhBWvfR+9dtmf/wAFTv2VPEH7L3/BtV+0ZdeOPso+JPxS8UWfj7xfHbTCa3s9Q1DxNpj/AGSNxwy28CwQbhwxiZhwRWh/wZU/8osvH3/ZVdR/9NGj18v/APBZn/g6B+AX/BRL/gmx8SPg54K8IfGDS/E/jD+zPsV1reladBYRfZtUs7yTzHhvpZBmO3cDbG2WKg4GSPL/APg3p/4OFvgv/wAEmf2L/E/w5+I3hj4oa1reteNbrxJBP4b06xubVLeWxsLdUZp7yFxIHtZCQEIwV+YkkBSd7JKySSS7JKySv0SSSLSd227tttvu27t/NtvTTsf03UV+QP8AxGrfssf9CD+0B/4I9I/+WdH/ABGrfssf9CD+0B/4I9I/+WdSM/X6ivyB/wCI1b9lj/oQf2gP/BHpH/yzo/4jVv2WP+hB/aA/8Eekf/LOgD9d9S1O20axkury4gtLaIZeWaQRogzjljwOTU9fg9+3r/wc7fsSf8FFv2XfEfwm8feAv2lYvD3iLyXa40zTNHgu7SaGVZYpY2OospKuoOHVlIyCDXxH+wR/wcofEf8A4Jq/EL/hE9B1/wAUfHH4BWcixadpXjm2j0zXNPtsDC280U90ICgwojMksOF+VIy3AB/V7RXiH/BPn9vrwV/wUk/Zw0z4m+BLLxRp+jag7QPb67pMthPDMoG9FZgYp1GceZA8kecjduDKPb6ACiiigAr8Af8Ag+c/5td/7mv/ANwtfv8AV+AP/B85/wA2u/8Ac1/+4WgD9f8A/gk7/wAosv2af+yVeF//AE0Wtd/8RPj18Pvh98dvh14G8R6tp9p48+If9p/8IhZTWryT6h9it1nvfKkCFY9kLKzbmXcDgZPFcB/wSd/5RZfs0/8AZKvC/wD6aLWvN/21P2bvG/xN/wCCuH7E/wAQ9C0C41Dwb8MP+E6/4SjU0liWPSft+jQ29puVmDt5kqso2K2COcDmgD7Hrz//AIRb4Zf8NT/235Hg/wD4XL/win2Hzt8H9v8A9gfbN+3bnzvsf2rnONnmd91egV8Mf8K91/8A4iXP+Er/ALD1j/hFv+GZf7J/tj7FJ/Z/2z/hKvN+zefjy/O8v5/L3btvOMc0Afc9ef8Aw78LfDLSfjt8RdT8LQeD4/iPq39mf8JxJpzwHVpfLt2XT/t4Q+YMQFhF5gGUzt4r0Cvhj9gr4e6/4c/4LV/t96/qOh6xYaF4k/4V5/ZGo3FlJFaap5GgzxzeRKwCS+W5CvsJ2sQDg0AfY/xZ8e+HvhV8LPE3ijxbd29h4U8N6VdaprNzPE0sVvZQQvLPI6KGLKsauSACSBgA9KPhP498PfFX4WeGvFHhK7t7/wAKeJNKtdU0a5giaKK4sp4UlgkRGClVaNkIBAIBwQOled/8FFPhvrfxk/4J9/HXwh4ZsJNV8SeK/h7r+j6VZRuqPeXdxptxDDEGYhQWkdVyxAGeSBR/wTr+G+t/Bv8A4J9/Arwh4msJNK8SeFPh7oGj6rZSOrvZ3dvptvDNEWUlSVkRlypIOOCRQB+Mv/B85/za7/3Nf/uFr9f/APgk7/yiy/Zp/wCyVeF//TRa1+QH/B85/wA2u/8Ac1/+4Wv1/wD+CTv/ACiy/Zp/7JV4X/8ATRa0AfP/APwTw/5Tr/8ABRX/ALpr/wCo9cV9/wBef/Dv49fD74g/Hb4i+BvDmrafd+PPh5/Zn/CX2UNq8c+n/bbdp7LzZCgWTfCrMu1m2gYODxXoFAHwB/wa4/8AKCj4Gf8Acf8A/Uh1Oj/gvp/zZX/2dV4G/wDb6vr/APZc+PXw+/ad+BOheOfhZq2n654D1z7R/Zd7Y2r20E3lXEsE22N0RlxNHKpyoyVJ5ByT4+/Hr4ffAf8A4Qr/AIWBq2n6V/wmfiux8LeG/tVq8/2zWrrzPssEe1G2SNskw7bVGDlhmgD0CvgD/ggX/wA3qf8AZ1Xjn/2xr7/rz/4BfHr4ffHj/hNf+Ff6tp+q/wDCGeK77wt4k+y2rwfY9atfL+1QSbkXfIu+PLruU5GGOKAPkD/g6O/5QUfHP/uAf+pDplff9ef/ALUfx6+H37MXwJ13xz8U9W0/Q/Aeh/Z/7Uvb61e5gh824igh3Rojs2ZpIlGFOCwPAGR6BQB8Af8ABPD/AJTr/wDBRX/umv8A6j1xR/wX0/5sr/7Oq8Df+31fX/w7+PXw++IPx2+Ivgbw5q2n3fjz4ef2Z/wl9lDavHPp/wBtt2nsvNkKBZN8Ksy7WbaBg4PFHx9+PXw++A//AAhX/CwNW0/Sv+Ez8V2Phbw39qtXn+2a1deZ9lgj2o2yRtkmHbaowcsM0AegV8Af8EC/+b1P+zqvHP8A7Y19/wBef/AL49fD748f8Jr/AMK/1bT9V/4QzxXfeFvEn2W1eD7HrVr5f2qCTci75F3x5ddynIwxxQB8gf8ABfT/AJsr/wCzqvA3/t9X3/Xn/wAffj18PvgP/wAIV/wsDVtP0r/hM/Fdj4W8N/arV5/tmtXXmfZYI9qNskbZJh22qMHLDNegUAfzA/8ABlT/AMpTfH3/AGSrUf8A076PX9P1fzA/8GVP/KU3x9/2SrUf/Tvo9fov8Qv+Dhf9pzwd4+1zSLH/AIJtfHjW7HStQns7fUYJdW8q/jjkZFmTGisNrgBhhmGG6nrQB+s9FfkD/wARHf7U/wD0jH/aA/7/AGr/APyjo/4iO/2p/wDpGP8AtAf9/tX/APlHQB+v1FfkD/xEd/tT/wDSMf8AaA/7/av/APKOj/iI7/an/wCkY/7QH/f7V/8A5R0Afr9RX5A/8RHf7U//AEjH/aA/7/av/wDKOj/iI7/an/6Rj/tAf9/tX/8AlHQB+v1RX9/BpVjPdXU0NtbW0bSzTSuEjiRRlmZjwAACSTwMV+Q3/ER3+1P/ANIx/wBoD/v9q/8A8o6g1L/g4p/ae1nTrizvP+CX/wAebu0u42hngmbVpI5kYEMjKdCwykEgg8EGgD9T/gP+0p8Pf2o/B7+IPhv428LeOtEima2kvdC1OG/hilXqjNGxCt3wcHBB6EV21fyAftI/Hj4nf8E8/wBo/TPi78F/gB8c/wBiSXW5mSbTdc1C9u9E1d1O/wAmKK90+3DxDLM0ErToMjasYAFfvd/wQK/4LCfE3/gqf8Jbm6+IHwZ1rwvJo8W3/hNdPjEXhrXpFIUpCsziVZueVi85BglnjyqEA/RGiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKAPgD/g6O/5QUfHP/uAf+pDpleAf8GVP/KLLx9/2VXUf/TRo9e//wDB0d/ygo+Of/cA/wDUh0yvAP8Agyp/5RZePv8Asquo/wDpo0egD9fqKKKACiiigDxT/goTL8cU/ZR8SL+znF4Yk+LExgh0p9edVtbdGlVZpQG+RpEjLMokymRyG+63wr+xH/wbTaVL8TU+MX7YPjO8/aM+MF4y3D2mpTST+H9NYciPy5ADdKhJCq6pAAcCDgNX6q0UAQaZpltomm29lZW8FpZ2kSwwQQxiOKGNQFVFUcKoAAAHAAqeiigAooooAK/AH/g+c/5td/7mv/3C1+/1fgD/AMHzn/Nrv/c1/wDuFoA/X/8A4JO/8osv2af+yVeF/wD00WtWPj/+3Ha/Af8Abd/Z9+C03h241K6+Pf8Awkf2fVUvBFHo39j2Ed426IoTL5ofYMMu0jPPSq//AASd/wCUWX7NP/ZKvC//AKaLWq/7SH7Dl18eP2+v2bPjTD4it9NtfgJ/wk/2jSnszLJrP9sadHZrtlDgReUU3nKtuBxx1oA+h68n/wCGyfCP/Dcv/DPnl6x/wnn/AAgv/Cw/M+zr/Z/9m/2h/Z+PN37vO87+DZjbzuzxXrFfLH/DG3i7/h9r/wANB+Zo/wDwgf8Awo//AIV55f2hv7Q/tL+3v7Qz5Wzb5Pk/x787uNuOaAPqevJ/hF+2T4R+Nf7UXxf+EWjx6wviv4Jf2N/wkL3FuqWj/wBq2j3dr5DhyXxGh35VdrYAz1r1ivlj9k39jbxd8FP+Cmv7W3xd1iTR28KfG3/hDv8AhHkt7hnu0/srSpbS689CgCZkcbMM25ck46UAe1/tTfHGL9mL9mP4jfEqfTpNYg+HnhfU/E0lhHMIXvVsrSW5MIchghcR7dxBxnOD0o/ZZ+OMX7Tv7Mfw5+JUGnSaPB8Q/C+meJo7CSYTPZLe2kVyIS4Chygk27gBnGcDpVf9rz4HS/tO/snfFD4awajHo8/xD8Jar4Zjv5ITMlk17Zy2wmKAqXCGTdtBGcYyOtH7IfwOl/Zi/ZO+F/w1n1GPWJ/h54S0rwzJfxwmFL1rKzitjMEJYoHMe7aScZxk9aAPxF/4PnP+bXf+5r/9wtfr/wD8Enf+UWX7NP8A2Srwv/6aLWvyA/4PnP8Am13/ALmv/wBwtfr/AP8ABJ3/AJRZfs0/9kq8L/8ApotaAPn/AP4J4f8AKdf/AIKK/wDdNf8A1Hrivv8Ar44/Yr/Zu8b/AAy/4K4ftsfEPXdAuNP8G/E//hBf+EX1N5Ymj1b7Bo01vd7VVi6+XKyqd6rknjI5r7HoA+AP+DXH/lBR8DP+4/8A+pDqdH/BfT/myv8A7Oq8Df8At9Xcf8EAP2bvG/7I3/BI/wCEvw8+I2gXHhfxl4e/tj+0NMnlilktvO1m/uIstGzId0UsbcMeG5wcij/gsP8As3eN/wBo7/hln/hCdAuNe/4QP9oXwn4w17ypYo/7O0m0+1/abpt7LuVPMTKrlju4U0AfY9fAH/BAv/m9T/s6rxz/AO2Nff8AXxx/wR4/Zu8b/s4/8NTf8JtoFxoP/CeftC+LPGGg+bLFJ/aOk3f2T7NdLsZtqv5b4VsMNvKigDh/+Do7/lBR8c/+4B/6kOmV9/18cf8ABf8A/Zu8b/tc/wDBI/4tfDz4c6BceKPGXiH+x/7P0yCWKKS58nWbC4lw0jKg2xRSNyw4XjJwK+x6APgD/gnh/wAp1/8Agor/AN01/wDUeuKP+C+n/Nlf/Z1Xgb/2+ruP2K/2bvG/wy/4K4ftsfEPXdAuNP8ABvxP/wCEF/4RfU3liaPVvsGjTW93tVWLr5crKp3quSeMjmj/AILD/s3eN/2jv+GWf+EJ0C417/hA/wBoXwn4w17ypYo/7O0m0+1/abpt7LuVPMTKrlju4U0AfY9fAH/BAv8A5vU/7Oq8c/8AtjX3/Xxx/wAEeP2bvG/7OP8Aw1N/wm2gXGg/8J5+0L4s8YaD5ssUn9o6Td/ZPs10uxm2q/lvhWww28qKAOH/AOC+n/Nlf/Z1Xgb/ANvq+/6+OP8AgsP+zd43/aO/4ZZ/4QnQLjXv+ED/AGhfCfjDXvKlij/s7SbT7X9pum3su5U8xMquWO7hTX2PQB/MD/wZU/8AKU3x9/2SrUf/AE76PX9P1fzA/wDBlT/ylN8ff9kq1H/076PX9P1ABRRRQAUUUUAFFFFABVLxJaX1/wCHb+DTLuPT9SmtpI7S6kh85LaUqQkhTI3hWwduRnGM1dooA/K39lf/AINndL1/43XHxe/bA+Iuo/tMfEmabfDaX/mR6DZqrEojQsczoOqxYjgUMV8lhg1+pGiaHZeGdHtdO02ztdP0+xiWC2tbaJYobeNRhURFACqAAAAMACrVFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFAHwB/wdHf8AKCj45/8AcA/9SHTK8A/4Mqf+UWXj7/squo/+mjR69/8A+Do7/lBR8c/+4B/6kOmV4B/wZU/8osvH3/ZVdR/9NGj0Afr9RRRQAUUUUAFFFFABRRRQAUUUUAFfgD/wfOf82u/9zX/7ha/f6vwB/wCD5z/m13/ua/8A3C0Afr//AMEnf+UWX7NP/ZKvC/8A6aLWuH/bD/aw8b/CL/gqn+xx8LtC1K3tvBvxg/4TX/hKLR7SKWS9/s3SIbq02yMpeLZK7E7CNwODkV3H/BJ3/lFl+zT/ANkq8L/+mi1rt/iZ+yf4I+Lv7Qvwx+KOu6bcXPjL4P8A9q/8IvdpdyxR2X9pWy2t3ujVgku+JFA3g7SMjBoA9Ir5o/4ba1//AIfGf8M4/wBk6P8A8It/wpr/AIWT/aeJP7Q+2f23/Z3kfe8vyfL+b7u7d/FjivpevF/+GJdA/wCHh/8Aw0d/a2sf8JT/AMK6/wCFbf2ZmP8As/7H/af9o+f93zPO8z5fvbdv8OeaAPaK+aP2Xv22tf8Ajp/wUY/ak+DWo6To9noXwJ/4RP8Asi9txJ9r1D+19MlvJvP3MU+R0CpsVflJzk819L14v8Ef2JdA+Bf7YHxx+MunatrF5rvx2/sH+17K4Mf2TT/7Isns4fI2qH+dHLPvZvmAxgcUAXP2/vixrfwF/YP+Nnjrwzcx2fiTwX4C13XdKuJIVmSC7tdOnnhcowKsBIinawIOMEYo/YB+LGt/Hr9g/wCCfjrxNcx3niTxp4C0LXdVuI4VhSe7utOgnmcIoCqDI7HaoAGcAYruPjJ8J9E+PXwh8VeBfE1tJeeG/Gmj3eharbxzNC89pdQvBMgdSGUmN2G5SCM5BzR8G/hPonwF+EPhXwL4ZtpLPw34L0e00LSreSZpngtLWFIIULsSzERoo3MSTjJOaAPwl/4PnP8Am13/ALmv/wBwtfr/AP8ABJ3/AJRZfs0/9kq8L/8Apota/ID/AIPnP+bXf+5r/wDcLX6//wDBJ3/lFl+zT/2Srwv/AOmi1oAsfAD9uO1+PH7bv7QXwWh8O3Gm3XwE/wCEc+0aq94JY9Z/tiwkvF2xBAYvKCbDlm3E546V7vXwB/wTw/5Tr/8ABRX/ALpr/wCo9cV9/wBAHhH/AATN/bjtf+CkX7EXgn402Xh248J2vjP7ds0qe8F3Ja/Zb+5szmUIgbcbcv8AdGA+OcZJ+3X+3Ha/sQ/8Kc+1eHbjxF/wt74oaJ8M4fKvBbf2ZJqXn7btso3mLH5JzGNpbd94Yr53/wCDXH/lBR8DP+4//wCpDqdH/BfT/myv/s6rwN/7fUAff9eEfsK/tx2v7b3/AAuP7L4duPDv/Cofihrfwzm828Fz/acmm+Ruu1wi+WsnnDEZ3Fdv3jmvd6+AP+CBf/N6n/Z1Xjn/ANsaAPoj/gpl+3Ha/wDBN39iLxt8ab3w7ceLLXwZ9h36VBeC0kuvtV/bWYxKUcLtNwH+6chMcZyPd6+AP+Do7/lBR8c/+4B/6kOmV9/0AeEfAD9uO1+PH7bv7QXwWh8O3Gm3XwE/4Rz7Rqr3glj1n+2LCS8XbEEBi8oJsOWbcTnjpR+3X+3Ha/sQ/wDCnPtXh248Rf8AC3vihonwzh8q8Ft/Zkmpeftu2yjeYsfknMY2lt33hivnf/gnh/ynX/4KK/8AdNf/AFHrij/gvp/zZX/2dV4G/wDb6gD7/rwj9hX9uO1/be/4XH9l8O3Hh3/hUPxQ1v4ZzebeC5/tOTTfI3Xa4RfLWTzhiM7iu37xzXu9fAH/AAQL/wCb1P8As6rxz/7Y0AfRH7df7cdr+xD/AMKc+1eHbjxF/wALe+KGifDOHyrwW39mSal5+27bKN5ix+ScxjaW3feGK93r4A/4L6f82V/9nVeBv/b6vv8AoA/mB/4Mqf8AlKb4+/7JVqP/AKd9Hr+n6v5gf+DKn/lKb4+/7JVqP/p30ev6fqACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooA+AP+Do7/lBR8c/+4B/6kOmV4B/wZU/8osvH3/ZVdR/9NGj17//AMHR3/KCj45/9wD/ANSHTK8A/wCDKn/lFl4+/wCyq6j/AOmjR6AP1+ooooAKKKKACiiigAooooAKKKKACvwB/wCD5z/m13/ua/8A3C1+/wBX4A/8Hzn/ADa7/wBzX/7haAP1/wD+CTv/ACiy/Zp/7JV4X/8ATRa14P8A8FBtWurP/guT/wAE87WG5uIrW8/4WR9ohSQrHPt8PwFdyjhsHkZ6GveP+CTv/KLL9mn/ALJV4X/9NFrXf/ET49fD74ffHb4deBvEerafaePPiH/af/CIWU1q8k+ofYrdZ73ypAhWPZCys25l3A4GTxQB6BXxR/wuvxd/xEZf8K5/4SPWP+ED/wCGcP8AhJP7A+0t/Z/9pf8ACT/Z/tnlZ2+d5P7vfjO3jpX2vXn/APwi3wy/4an/ALb8jwf/AMLl/wCEU+w+dvg/t/8AsD7Zv27c+d9j+1c5xs8zvuoA9Ar4o/Yb+Nfi7xz/AMFjv26vBuseI9Y1Pwp4F/4QH/hHtJuLlpLTRfteiTTXXkRk4j82RQ74+8wBNfa9ef8Aw78LfDLSfjt8RdT8LQeD4/iPq39mf8JxJpzwHVpfLt2XT/t4Q+YMQFhF5gGUzt4oA4T/AIKr3cth/wAEvP2kZ4JZIJ4Phb4nkjkjYq8bDSbohgRyCDzkUf8ABKi7lv8A/gl5+zdPPLJPPP8AC3wxJJJIxZ5GOk2pLEnkknnJr1v4s+PfD3wq+FnibxR4tu7ew8KeG9KutU1m5niaWK3soIXlnkdFDFlWNXJABJAwAelHwn8e+Hvir8LPDXijwld29/4U8SaVa6po1zBE0UVxZTwpLBIiMFKq0bIQCAQDggdKAPwh/wCD5z/m13/ua/8A3C1+v/8AwSd/5RZfs0/9kq8L/wDpota/ID/g+c/5td/7mv8A9wtfr/8A8Enf+UWX7NP/AGSrwv8A+mi1oAr/ALN/7Dl18B/2+v2k/jTN4it9Stfj3/wjH2fSkszFJo39j6dJZtulLkS+aX3jCrtAxz1r6Hr5I/Y8/aw8b/F3/gqn+2P8Ltd1K3ufBvwf/wCEK/4Re0S0iiksv7S0ia6u90iqHl3yopG8naBgYFfW9AHzx/wSk/Ycuv8Agm7+wL4C+C174it/Fl14M/tDfqsFmbSO6+1ajdXgxEXcrtFwE+8clM8ZwD/goV+w5dftvf8ACjvsviK38O/8Kh+LWg/EybzbM3P9px6b9o3Wi4dfLaTzhiQ7gu37pzXH/wDBCn9rDxv+3F/wSs+FvxR+I2pW+r+MvFH9rf2hdwWkVpHL9n1e9tYsRxqqLiKGMcAZIyeSaP8AgrR+1h43/ZS/4Zn/AOEJ1K307/hZPx68LeBNe820iuPtWk332r7TCu9T5bN5SYdcMuOCM0AfW9fPH/BPX9hy6/Yh/wCF4/avEVv4i/4W98Wte+JkPlWZtv7Mj1L7PttGy7eY0fknMg2ht33Rivoevkj/AIJL/tYeN/2rf+GmP+E21K31H/hW3x68U+BNB8q0it/suk2P2X7NC2xR5jL5r5dss2eScUAdh/wVb/Ycuv8AgpF+wL49+C1l4it/Cd14z/s/Zqs9mbuO1+y6ja3hzEHQtuFuU+8MF884wfoevkj/AILrftYeN/2Hf+CVnxT+KPw51K30jxl4X/sn+z7ue0iu44vtGr2VrLmORWRsxTSDkHBORyBX1vQB88fs3/sOXXwH/b6/aT+NM3iK31K1+Pf/AAjH2fSkszFJo39j6dJZtulLkS+aX3jCrtAxz1o/4KFfsOXX7b3/AAo77L4it/Dv/Cofi1oPxMm82zNz/acem/aN1ouHXy2k84YkO4Lt+6c1x/7Hn7WHjf4u/wDBVP8AbH+F2u6lb3Pg34P/APCFf8IvaJaRRSWX9paRNdXe6RVDy75UUjeTtAwMCj/grR+1h43/AGUv+GZ/+EJ1K307/hZPx68LeBNe820iuPtWk332r7TCu9T5bN5SYdcMuOCM0AfW9fPH/BPX9hy6/Yh/4Xj9q8RW/iL/AIW98Wte+JkPlWZtv7Mj1L7PttGy7eY0fknMg2ht33Rivoevkj/gkv8AtYeN/wBq3/hpj/hNtSt9R/4Vt8evFPgTQfKtIrf7LpNj9l+zQtsUeYy+a+XbLNnknFAHYf8ABQr9hy6/be/4Ud9l8RW/h3/hUPxa0H4mTebZm5/tOPTftG60XDr5bSecMSHcF2/dOa+h6+SP+CtH7WHjf9lL/hmf/hCdSt9O/wCFk/Hrwt4E17zbSK4+1aTffavtMK71Pls3lJh1wy44IzX1vQB/MD/wZU/8pTfH3/ZKtR/9O+j1/T9X8wP/AAZU/wDKU3x9/wBkq1H/ANO+j1/T9QAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFABRRRQAUUUUAFFFFAHwB/wAHR3/KCj45/wDcA/8AUh0yvAP+DKn/AJRZePv+yq6j/wCmjR69/wD+Do7/AJQUfHP/ALgH/qQ6ZXgH/BlT/wAosvH3/ZVdR/8ATRo9AH6/UUUUAFFFFABRRRQAUUUUAFFFFABX4A/8Hzn/ADa7/wBzX/7ha/f6vwB/4PnP+bXf+5r/APcLQB+v/wDwSd/5RZfs0/8AZKvC/wD6aLWvN/21P2bvG/xN/wCCuH7E/wAQ9C0C41Dwb8MP+E6/4SjU0liWPSft+jQ29puVmDt5kqso2K2COcDmvSP+CTv/ACiy/Zp/7JV4X/8ATRa1Y+P/AO3Ha/Af9t39n34LTeHbjUrr49/8JH9n1VLwRR6N/Y9hHeNuiKEy+aH2DDLtIzz0oA93r4Y/4V7r/wDxEuf8JX/Yesf8It/wzL/ZP9sfYpP7P+2f8JV5v2bz8eX53l/P5e7dt5xjmvuevJ/+GyfCP/Dcv/DPnl6x/wAJ5/wgv/Cw/M+zr/Z/9m/2h/Z+PN37vO87+DZjbzuzxQB6xXwx+wV8Pdf8Of8ABav9vvX9R0PWLDQvEn/CvP7I1G4spIrTVPI0GeObyJWASXy3IV9hO1iAcGvuevJ/hF+2T4R+Nf7UXxf+EWjx6wviv4Jf2N/wkL3FuqWj/wBq2j3dr5DhyXxGh35VdrYAz1oAp/8ABRT4b638ZP8Agn38dfCHhmwk1XxJ4r+Huv6PpVlG6o95d3Gm3EMMQZiFBaR1XLEAZ5IFH/BOv4b638G/+CffwK8IeJrCTSvEnhT4e6Bo+q2Ujq72d3b6bbwzRFlJUlZEZcqSDjgkV1H7U3xxi/Zi/Zj+I3xKn06TWIPh54X1PxNJYRzCF71bK0luTCHIYIXEe3cQcZzg9KP2WfjjF+07+zH8OfiVBp0mjwfEPwvpniaOwkmEz2S3tpFciEuAocoJNu4AZxnA6UAfh7/wfOf82u/9zX/7ha/X/wD4JO/8osv2af8AslXhf/00WtfkB/wfOf8ANrv/AHNf/uFr9f8A/gk7/wAosv2af+yVeF//AE0WtAHz/wD8E8P+U6//AAUV/wC6a/8AqPXFff8AXm/wz/ZP8EfCL9oX4nfFHQtNuLbxl8YP7K/4Si7e7lljvf7NtmtbTbGzFItkTsDsA3E5OTXpFAHwB/wa4/8AKCj4Gf8Acf8A/Uh1Oj/gvp/zZX/2dV4G/wDb6vrf9kf9k/wR+w7+z14f+F3w50240jwb4X+0/wBn2k93LdyRfaLmW6lzJIzO2ZZpDyTgHA4Ao/aN/ZP8EftW/wDCB/8ACbabcaj/AMK28Yaf470HyruW3+y6tY+Z9mmbYw8xV818o2VbPIOKAPSK+AP+CBf/ADep/wBnVeOf/bGvv+vN/wBnL9k/wR+yl/wnn/CE6bcad/wsnxhqHjvXvNu5bj7Vq195f2mZd7Hy1bykwi4VccAZoA+SP+Do7/lBR8c/+4B/6kOmV9/15v8Atcfsn+CP24v2evEHwu+I2m3Gr+DfFH2b+0LSC7ltJJfs9zFdRYkjZXXEsMZ4IyBg8E16RQB8Af8ABPD/AJTr/wDBRX/umv8A6j1xR/wX0/5sr/7Oq8Df+31fW/wz/ZP8EfCL9oX4nfFHQtNuLbxl8YP7K/4Si7e7lljvf7NtmtbTbGzFItkTsDsA3E5OTR+0b+yf4I/at/4QP/hNtNuNR/4Vt4w0/wAd6D5V3Lb/AGXVrHzPs0zbGHmKvmvlGyrZ5BxQB6RXwB/wQL/5vU/7Oq8c/wDtjX3/AF5v+zl+yf4I/ZS/4Tz/AIQnTbjTv+Fk+MNQ8d695t3LcfatWvvL+0zLvY+WreUmEXCrjgDNAHyR/wAF9P8Amyv/ALOq8Df+31ff9eb/ALRv7J/gj9q3/hA/+E20241H/hW3jDT/AB3oPlXctv8AZdWsfM+zTNsYeYq+a+UbKtnkHFekUAfzA/8ABlT/AMpTfH3/AGSrUf8A076PX9P1fzA/8GVP/KU3x9/2SrUf/Tvo9f0/UAFFFFABRRRQB8t/t9/tpeLPhX8V/hj8E/hFZ6DqPxk+L11K9rcayjz6b4U0a1Ae+1a5hjeN5tinZFCJI/Mkb7+EKtwPjr4wfFD9gf8AbT/Z78HeJfiv4n+MvhL463194Z1BPEOh6RaXuiahb2huYbu0fTLO1HkuwKSRTrKVXayuNrbs7wrp8cf/AAch+Lb3xEsVvcyfA6xtfCBnKqbyEarK995OeWZHMe8DkKwJGDmqL7P+CjX/AAWh8JeJPDsq6r8Jv2SNO1O3udYij32Wp+MNQQW8llDJjZK1nbKGkZCfKlcIcMeFhVf2Ul9pzcuyhCU4NW6JqNk9/aTWvwJGIaSqx2UVFR7uc4Rkn52lJNrbkhJ21lzfoDRXzD8Cv209bvdb/aa1jx7e6FceAPgp4qk0fTb7w/4Z1D7b9nh0+3vLrz4UmupLmSE3KxF7eJAxgkbyxu2o74X/APBRbwv4J/Z6+FWvfGLxn4RtfEXxNtrO7trnwzo+rf2MsWoTBbCSQzxGayhkEkMYmvhArS7hhT8ikWmk11UGv+4ivH5tJ272dhyVm0+jkv8AwF2l92l+11ex9OUV5f8AA/8AbL+HH7RTeLh4V1+eZvAskaa4NS0q80g2SSxedFOBeRReZbyRAyR3Ee6GRBuV2HNct4B/4Ka/BP4m32qw6R4vuni0rS31tby68P6nZWOr2Kzrb/atOuZrdIdSiaaSONHsnmDtLGFLb0yNpOz7X+Vm7+lk3fsmxLVcy2vb53tb1u0vXQ95orwj4Yf8FLvgz8ZvEnhHR/DXibVdT1LxtdXdhp0A8MatEYLq1NyJ7a8L2yiwnX7HdYiuzC7CCQqrBTXuGq2/2zS7mLyLe682Jk8mc4imyCNrnDfKeh4PB6HpTneMW7d/wCNnLluT0V8Af8O5P+rBf2AP/Ch//BCvZf2KP2U/+FEfEXU9T/4Zn/Zk+C32rTzbf2t8N9U+1ajeZkRvs8q/2Fp+ITt3E+a/zIvyfxLUVd2ZMnZaH0fr2vWXhbQr3U9SuoLHTtOge6urmZwkdvEilndieAqqCSewFfE3"
                     +
                    "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",
                fileName=
                    "U:/Documents/PHD/Articles/1st paper/gfx/tank_wall.jpg")}));
      end HeatTransferTank_PCM;

      package PCM_NEW_with_Nat_Conv_Coupling_Area
        model PCM_inner

          import SI = Modelica.SIunits;
        type Specific_Heat_Capacity_Coefficient = Real(quantity="Specific_Heat_Capacity_Coefficien", unit="J/(kg.K2)");

        replaceable package Medium =
            ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
            Modelica.Media.Interfaces.PartialMedium;
            Medium.ThermodynamicState medium;

        /********************************************************************************************************************************/
        /**************** PCM PROPERTIES [1]  ****************/
        /********************************************************************************************************************************/

        // PARAMETERS
        parameter SI.KinematicViscosity nu = 6.5e-6;
        parameter SI.LinearExpansionCoefficient beta1 = 8e-3;
        parameter SI.SpecificHeatCapacity  c_s_PCM = G2*2905;  // in the paper 2905"J/(kg K)"
        parameter SI.SpecificHeatCapacity c_l_PCM = G2*2905; // in the paper 2905"J/(kg K)" same value for solid and liquid specific heat capacity
        parameter SI.SpecificEnthalpy lambda = G3*224360 "J/kg"; // latent heat
        parameter SI.Density rho_s_PCM = G1*785.3;
        parameter SI.Density rho_l_PCM = G1*769.2;
        parameter SI.Density rho_PCM_m = rho_s_PCM + (rho_l_PCM-rho_s_PCM)*(T_PCM_m - T_m1)/(T_m2 -T_m1)
            "linearly interpolated density at T = T_m";

        parameter SI.ThermalConductivity k_s_PCM = Go*0.24 "W/(m K)";
        parameter SI.ThermalConductivity k_l_PCM = Go*0.18 "W/(m K)";

        parameter SI.DynamicViscosity mu = nu*rho_l_PCM;
        parameter Real Pr= mu*c_l_PCM/k_l_PCM; // Prandtl number
        parameter SI.SpecificHeatCapacity  c_PCM_m=(2*lambda-2*(c_l_PCM-c_s_PCM)*(T_m2-T_PCM_m) + c_s_PCM*(T_PCM_m-T_m1) + c_l_PCM*(T_m2-T_PCM_m))/(T_m2-T_m1);
        parameter SI.SpecificHeatCapacity a1= c_s_PCM; // coefficient in c_PCM calculation
        parameter SI.SpecificHeatCapacity  c = c_PCM_m; // coefficient in c_PCM calculation
        parameter SI.SpecificHeatCapacity c_PCM_m1=2*lambda/(T_m2-T_m1) + c_s_PCM;
        parameter SI.Acceleration g=9.81;

        // VARIABLES
        SI.SpecificHeatCapacity c_PCM;
        SI.Density rho_PCM;
        SI.ThermalConductivity k_PCM "W/(m K)";  // effective conductivity of the PCM
        SI.ThermalConductivity k_PCM_m
            "linearly interpolated thermal conductivity at T = T_m";
        Real k_eff_ratio "W/(m K)";
        Real k_eff_ratio2 "W/(m K)";

        /*******************************************************************************************************************************/
        /**************** PCM TEMPERATURE [1] ****************/
        /*******************************************************************************************************************************/

        // PARAMETERS
        parameter SI.Temperature T_PCM_start = 293 "ambient temperature";
        parameter SI.Temperature T_PCM_m = 328.15 "55 C + 273"; // melting temperature
        /* Definition of an assymmetrical T-melting range (i.e. assymmetrical triangle for c_PCM_eff during phase change */
        parameter SI.Temperature T_m1 = T_PCM_m - 3; // low value of the temperature melting range
        parameter SI.Temperature T_m2 = T_PCM_m + 1; // high value of the temperature melting range

        // VARIABLES
        SI.Temperature T_PCM(start=T_PCM_start, fixed=true);
        Specific_Heat_Capacity_Coefficient b_PCM; // T-dependent coefficient in c_PCM calculation
        Specific_Heat_Capacity_Coefficient d_PCM; // T-dependent coefficient in c_PCM calculation
        SI.ThermalResistance R_th "K/W";
        SI.Volume V_dx "PCM discretized volumes in [m3]";
        SI.Mass M_dx; // PCM mass in the discretized control volume (varies with density up to 7%)
        //SI.SpecificEnthalpy h_PCM;
        Real Ra1; // Rayleigh number
        Real Ra2; // Rayleigh number
        // Real Ra33; // Rayleigh number

        /*******************************************************************************************************************************/
        /**************** GEOMETRY ****************/
        /*******************************************************************************************************************************/

        // OUTER PARAMETERS
        //outer parameter SI.Length D "tank diameter in [m]";

        outer SI.Length d "tank diameter in [m]";
        outer SI.Length  L "tank length in [m]";
        outer SI.Length  Lc "length scale in the Rayleigh number [2]";
        outer SI.CoefficientOfHeatTransfer  h_i;
        outer SI.Length dx_PCM;
        /*
 parameter SI.Length d=0.34 "tank diameter in [m]";
 parameter SI.Length  L=0.8 "tank length in [m]";
 parameter SI.Length Lc=0.11 "length scale in the Rayleigh number [2]";
 parameter SI.CoefficientOfHeatTransfer  h_i=150;
*/
        //outer parameter SI.Length t_PCM "total PCM-layer thickness";
         SI.Length dx_inner= dx_PCM/2
            "discretization step for the inner surface = 0.1 mm";
        SI.Area A_hex "total area of heat transfer for a cylinder";

        /*******************************************************************************************************************************/
        /**************** BOOLEAN ****************/
        /*******************************************************************************************************************************/
        // VARIABLES
        Boolean SOLID;
        Boolean TWO_PHASE1;
        Boolean TWO_PHASE2;
        Boolean LIQUID;
        outer parameter Boolean Convection; // if false only conduction is considered in the PCM

         input SI.Temperature Delta_T1;
        //  input SI.Temperature Delta_T11;

        // input SI.Temperature Delta_T22;
         input Real Delta_T22;
         // input Real Delta_T33;

        SI.Temperature delta_T;
        //SI.Temperature Delta_T33;

        Real Ra11;
        //Real Ra22;
        //Real Ra33;
        //Real k_eff_ratio33;
        Real k_eff_ratio11;
        // Real k_eff_ratio33 "W/(m K)";
        //Real k_eff_ratio22;
        Real k_Test;
        Real k_Test1;

          Ports.HeatFlow2 A
                annotation (Placement(transformation(extent={{-108,-12},{-82,14}})));

          Ports.HeatFlow  B annotation (Placement(transformation(extent={{90,-10},{110,10}}),
                iconTransformation(extent={{80,-24},{120,20}})));

        /********************************************************************************************************************************/
        // EXTRA PARAMETERS AND VARIABLES FROM "InnerWallCell"
        /********************************************************************************************************************************/

          SI.Pressure p "Pressure of the hydrogen";
          SI.CoefficientOfHeatTransfer h_H2 "Heat transfer coefficient";
          Real Tau "Dimenasionless time";
          Real Ra "Dimensionless hydrogen properties number";
          Real beta "Thermal expansion coefficient";
          Real v "kinematic viscosity";
          Real a "thermal diffusivity of hydrogen";
          Real Nu "Dimensionless heat transfer number";
        //  outer Real   y1;
        SI.Length Radius_Local;
        // SI.Temperature   delta_T ;

          /****************** Tables *******************/
        /*  Modelica.Blocks.Tables.CombiTable1Ds tank_prop(
    tableName="properties",
    tableOnFile=true,
    table=[1,850,15,481; 2,2700,236,900; 3,1286,1.17,1578; 4,1374,1.14,1075],
    fileName=
        "C:/Users/edro/Documents/Dymola/External files/Lookuptables/HeatTransferProperties.txt")
    annotation (Placement(transformation(extent={{-80,0},{-60,20}})));
*/
        //      outer parameter SI.Temperature  T_amb "Ambient temperature";
        // SI.Temperature T "hydrogen temperature";

        //initial equation
        //T=T_amb;

        //   h_PCM = c_PCM*(T_PCM-273.15);

        outer parameter Real Go "multiplying factor for k_s and k_l";
        outer parameter Real G1 "multiplying factor for rho_s and rho_l";
        outer parameter Real G2 "multiplying factor for cp_s and cp_l";
        outer parameter Real G3 "multiplying factor for lambda";

        /************************************************************************************/
        outer Real G_A_inner
            "multiplying factor for inner heat transfer area at the H2/PCM interface";
        outer Real G_A_PCM
            "multiplying factor for inner heat transfer area in the R_th formule, so PCM_inner.B.Q is equal to PCM[1].A.Q as it should be";
        /************************************************************************************/

        equation
        (A_hex*dx_inner)*rho_PCM*c_PCM*der(T_PCM) = A.Q + B.Q; // energy balance in [W]

        //h_PCM = f*lambda + c_PCM*PCM;

        A_hex = (2*Modelica.Constants.pi*d/2*L + 2*Modelica.Constants.pi*(d/2)^2)
            "total area of heat transfer for a cylinder";

        Radius_Local = d/2;

        A.Q = (A.T-T_PCM)/(1/(h_H2*A_hex*G_A_inner));
        B.Q = G_A_PCM*(B.T-T_PCM)/(R_th)
            "not R_th/2 because dx_inner = dx_PCM/2 already";
        R_th = dx_inner/(A_hex*k_PCM); // thermal resistance

        b_PCM= (c_PCM_m-c_s_PCM)/(T_PCM_m-T_m1); // T-dependent coefficient in c_PCM calculation
        d_PCM= (c_PCM_m - c_l_PCM)/(T_m2-T_PCM_m); // T-dependent coefficient in c_PCM calculation

        M_dx = rho_PCM*A_hex*dx_inner;
        V_dx = A_hex*dx_inner; // discretized volume

        if Convection == true then
           k_eff_ratio = max(0.386*(Pr/(0.861 + Pr))^(1/4)*Ra1^(1/4), 1);
        k_PCM_m = k_s_PCM + (k_l_PCM*k_eff_ratio2-k_s_PCM)*(T_PCM_m - T_m1)/(T_m2 -T_m1); // linearly interpolated thermal conductivity at T = T_m
         k_eff_ratio2 = max(0.386*(Pr/(0.861 + Pr))^(1/4)*Ra2^(1/4),1);

        else
          k_eff_ratio = 1;
           k_eff_ratio2 = 1;

        k_PCM_m = k_s_PCM + (k_l_PCM-k_s_PCM)*(T_PCM_m - T_m1)/(T_m2 -T_m1);

        end if;

        // The specific heat capacity of the PCM is calculated depending on the physical condition (i.e. solid, two-phase: melting, liquid)
        if T_PCM < T_m1 then
          c_PCM = c_s_PCM;
          rho_PCM = rho_s_PCM;
          k_PCM = k_s_PCM;
          //Ra = 0;
          SOLID = true;
          TWO_PHASE1 = false;
          TWO_PHASE2 = false;
          LIQUID = false;

        else
          if T_m1<= T_PCM and T_PCM <= T_PCM_m then
          c_PCM= a1 + b_PCM*(T_PCM-T_m1);
          rho_PCM = rho_s_PCM + (rho_PCM_m - rho_s_PCM)*(T_PCM - T_m1)/(T_PCM_m-T_m1);
          k_PCM = k_s_PCM + (k_PCM_m -k_s_PCM)*(T_PCM - T_m1)/(T_PCM_m-T_m1);
          //k_PCM_m
          //Ra = 0;
          SOLID = false;
          TWO_PHASE1 = true;
          TWO_PHASE2 = false;
          LIQUID = false;

        //  k_eff_ratio = 1;

        //f=h_PCM/lambda;

        else if T_PCM_m < T_PCM and T_PCM< T_m2 then
          c_PCM = c - d_PCM*(T_PCM-T_PCM_m);
          rho_PCM = rho_PCM_m + (rho_l_PCM-rho_PCM_m)*(T_PCM - T_PCM_m)/(T_m2-T_PCM_m);
          k_PCM = k_PCM_m + ((k_l_PCM*k_eff_ratio)-k_PCM_m)*(T_PCM - T_PCM_m)/(T_m2-T_PCM_m);

          //Ra = 0;
        //  k_eff_ratio = 1;

          SOLID = false;
          TWO_PHASE1 = false;
          TWO_PHASE2 = true;
          LIQUID = false;

        //f=h_PCM/lambda;

        else
          c_PCM = c_l_PCM;
          rho_PCM = rho_l_PCM;

        k_PCM = k_l_PCM*k_eff_ratio;

        //f=1;

          SOLID = false;
          TWO_PHASE1 = false;
          TWO_PHASE2 = false;
          LIQUID = true;

        end if;
        end if;
        end if;

        // k_eff_ratio = 0.386*(Pr/(0.861 + Pr))^(1/4)*Ra1^(1/4);
         Ra1 =(g*beta1*abs(delta_T)*Lc^3)/(nu*k_l_PCM/(rho_l_PCM*c_l_PCM));

        Ra11 =(g*beta1*abs(Delta_T1)*Lc^3)/(nu*k_l_PCM/(rho_l_PCM*c_l_PCM));
         k_eff_ratio11 = max(0.386*(Pr/(0.861 + Pr))^(1/4)*Ra11^(1/4), 1);

        /*
 Ra22 =(g*beta1*abs(Delta_T22)*Lc^3)/(nu*k_l_PCM/(rho_l_PCM*c_l_PCM));
 k_eff_ratio22 = 0.386*(Pr/(0.861 + Pr))^(1/4)*Ra22^(1/4);
*/

        delta_T = Delta_T1;
        k_Test = k_l_PCM*k_eff_ratio2;
        k_Test1 = k_l_PCM*k_eff_ratio;

         // Ra2 =(g*beta1*(T_m2-T_PCM_start)*Lc^3)/(nu*k_l_PCM/(rho_l_PCM*c_l_PCM));
         Ra2 =max((g*beta1*(Delta_T22)*Lc^3)/(nu*k_l_PCM/(rho_l_PCM*c_l_PCM)), 0);
         /* Ra33 =max((g*beta1*(Delta_T33)*Lc^3)/(nu*k_l_PCM/(rho_l_PCM*c_l_PCM)),0); // to be compared with the correct one (i.e. Ra2)
 k_eff_ratio33 = max(0.386*(Pr/(0.861 + Pr))^(1/4)*Ra33^(1/4), 1); // to be compared with the correct one (i.e. k_eff_ratio_2)
*/
        // Counter is an Integer variable that keeps track of the actual position of the PCM block to calculate the correct heat transfer area
        // It is initialized with a value of 0 only for the "A" port of "PCM_inner"

        A.P=p;
        A.P=B.P;

        medium=Medium.setState_pT(p, A.T);
        Medium.thermalConductivity(medium)=a;
        Medium.dynamicViscosity(medium)/Medium.density(medium)=v;

        //calculation of rayleighs number taken from 'Natural convection cooling of
        //rectangular and cylindrical containers' by Wenxian Lin, S.W. Armfield
        Ra=g*beta*d^3*Medium.specificHeatCapacityCp(medium)*Medium.density(medium)^2*abs(A.T-T_PCM)/(v*k_PCM);
        beta=1/A.T;
        Tau=time/(d^2/(a*Ra^(1/2)));
        Nu=0.104*Ra^(0.352);
        /*
if A.m_flow < 0.00 then
  h=Nu*k_PCM/d;
elseif A.m_flow > 0.00 then
  h=h_i;
else
  h=50;
end if;
*/

         if A.m_flow >= 1e-15 then
         h_H2 = 150; // equal to "h_charging" in "HeatTransferTank"
         else

           h_H2 = 50;

         end if;

        // h_H2 = 100;

        A.Counter + 1 = B.Counter;

        /*
[1] Farid. M.M. et al. "Melting and solidifcation in multi.dimensional geometry and presence of more than one interface".
    Energy and Conversion Management, Vol. 38, No. 8 pp. 809-818, 1998. DOI: 10.1016/S0196-8904(97)00038-1.
    http://findit.dtu.dk/en/catalog/13046605 

[2] Incropera F. et al. "Introduction to Heat transfer". Wiley&Sons 6th Ed., pp. 592.593. 
    ISBN 13 978-0470-50196-2.
*/

          annotation (Diagram(graphics), Icon(graphics={                                Text(
                  extent={{-44,80},{44,54}},
                  lineColor={0,0,0},
                  textString="Inner"), Bitmap(extent={{-78,-64},{84,58}},
                    fileName=
                      "modelica://HySDeP/Graphics/PCM - rubrick cube figure for Dymola.jpg")}));
        end PCM_inner;

        model PCM

          import SI = Modelica.SIunits;
        type Specific_Heat_Capacity_Coefficient = Real(quantity="Specific_Heat_Capacity_Coefficien", unit="J/(kg.K2)");

        /********************************************************************************************************************************/
        /**************** PCM PROPERTIES [1]  ****************/
        /********************************************************************************************************************************/

        // PARAMETERS
        parameter SI.KinematicViscosity nu = 6.5e-6;
        parameter SI.LinearExpansionCoefficient beta = 8e-3;
        parameter SI.SpecificHeatCapacity  c_s_PCM = G2*2905;  // in the paper 2905"J/(kg K)"
        parameter SI.SpecificHeatCapacity c_l_PCM = G2*2905; // in the paper 2905"J/(kg K)" same value for solid and liquid specific heat capacity
        parameter SI.SpecificEnthalpy lambda = G3*224360 "J/kg"; // latent heat
        parameter SI.Density rho_s_PCM = G1*785.3;
        parameter SI.Density rho_l_PCM = G1*769.2;
        parameter SI.Density rho_PCM_m = rho_s_PCM + (rho_l_PCM-rho_s_PCM)*(T_PCM_m - T_m1)/(T_m2 -T_m1)
            "linearly interpolated density at T = T_m";
        parameter SI.ThermalConductivity k_s_PCM = Go*0.24 "W/(m K)";
        parameter SI.ThermalConductivity k_l_PCM = Go*0.18 "W/(m K)";
        parameter SI.DynamicViscosity mu = nu*rho_l_PCM;
        parameter Real Pr= mu*c_l_PCM/k_l_PCM; // Prandtl number
        parameter SI.SpecificHeatCapacity  c_PCM_m=(2*lambda-2*(c_l_PCM-c_s_PCM)*(T_m2-T_PCM_m) + c_s_PCM*(T_PCM_m-T_m1) + c_l_PCM*(T_m2-T_PCM_m))/(T_m2-T_m1);
        parameter SI.SpecificHeatCapacity a= c_s_PCM; // coefficient in c_PCM calculation
        parameter SI.SpecificHeatCapacity  c = c_PCM_m; // coefficient in c_PCM calculation
        parameter SI.SpecificHeatCapacity c_PCM_m1=2*lambda/(T_m2-T_m1) + c_s_PCM;
        parameter SI.Acceleration g=9.81;

        // VARIABLES
        SI.SpecificHeatCapacity c_PCM;
        SI.Density rho_PCM;
        SI.ThermalConductivity k_PCM "W/(m K)";  // effective conductivity of the PCM
        Real k_eff_ratio "W/(m K)";
        // SI.ThermalConductivity k_eff "W/(m K)";
        SI.ThermalConductivity k_PCM_m;
        Real k_eff_ratio2 "W/(m K)";

        /*******************************************************************************************************************************/
        /**************** PCM TEMPERATURE [1] ****************/
        /*******************************************************************************************************************************/

        // PARAMETERS
        parameter SI.Temperature T_PCM_start = 293 "55 C + 273"; // melting temperature
        parameter SI.Temperature T_PCM_m = 328.15 "55 C + 273"; // melting temperature
        /* Definition of an assymmetrical T-melting range (i.e. assymmetrical triangle for c_PCM_eff during phase change */
        parameter SI.Temperature T_m1 = T_PCM_m - 3; // low value of the temperature melting range
        parameter SI.Temperature T_m2 = T_PCM_m + 1; // high value of the temperature melting range

        // VARIABLES
        SI.Temperature T_PCM(start=T_PCM_start, fixed=true);
        Specific_Heat_Capacity_Coefficient b_PCM; // T-dependent coefficient in c_PCM calculation
        Specific_Heat_Capacity_Coefficient d_PCM; // T-dependent coefficient in c_PCM calculation
        SI.ThermalResistance R_th "K/W";
        SI.Volume V_dx "PCM discretized volumes in [m3]";
        SI.Mass M_dx; // PCM mass in the discretized control volume (varies with density up to 7%)
        //SI.SpecificEnthalpy h_PCM;
        Real Ra; // Rayleigh number
        Real Ra2; // Rayleigh number
        // Real Ra33; // Rayleigh number

        //SI.MassFraction f;

        /*******************************************************************************************************************************/
        /**************** GEOMETRY ****************/
        /*******************************************************************************************************************************/

        // OUTER PARAMETERS
        outer parameter SI.Length d "tank diameter in [m]";
        outer parameter SI.Length L "tank length in [m]";
        outer parameter SI.Length Lc "length scale in the Rayleigh number [2]";
        outer parameter SI.Length t_PCM "total PCM-layer thickness";
        outer parameter SI.Length dx_PCM "PCM discretization step in [m]";
        SI.Area A_hex "total area of heat transfer for a cylinder";

        SI.Length dx_TEST_PCM = dx_PCM;

        /*******************************************************************************************************************************/
        /**************** BOOLEAN ****************/
        /*******************************************************************************************************************************/
        // VARIABLES
        Boolean SOLID;
        Boolean TWO_PHASE1;
        Boolean TWO_PHASE2;
        Boolean LIQUID;
        outer parameter Boolean Convection;  // if false only conduction is considered

        SI.Length Radius_Local "Local radius";

        //input SI.Temperature delta_T2;
        SI.Temperature delta_T2;
        input SI.Temperature Delta_T1;
        // input SI.Temperature Delta_T22;
         input Real Delta_T22;
         //input Real Delta_T33;

        Real Ra11;
        Real Ra22;
        Real k_eff_ratio11;
        Real k_eff_ratio22;
        //Real k_eff_ratio33 "W/(m K)";
        Real k_Test;
        Real k_Test1;

          Ports.HeatFlow A
            annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
          Ports.HeatFlow B
            annotation (Placement(transformation(extent={{90,-10},{110,10}})));

        outer parameter Real Go "multiplying factor for k_s and k_l";
        outer parameter Real G1 "multiplying factor for rho_s and rho_l";
        outer parameter Real G2 "multiplying factor for cp_s and cp_l";
        outer parameter Real G3 "multiplying factor for lambda";
        outer parameter Real G_A_PCM
            "multiplying factor for the heat transfer area for each PCM volume";

        equation
        (A_hex*dx_PCM)*rho_PCM*c_PCM*der(T_PCM) = A.Q + B.Q; // energy balance in [W]

        //h_PCM = f*lambda + c_PCM*T_PCM;

        A_hex= (2*Modelica.Constants.pi*L*(d/2+ A.Counter*dx_PCM) + 2*Modelica.Constants.pi*(d/2+ A.Counter*dx_PCM)^2);

        Radius_Local = d/2+ A.Counter*dx_PCM;

        A.Q = (A.T-T_PCM)/(R_th/2);
        B.Q = (B.T-T_PCM)/(R_th/2);
        R_th = dx_PCM/(G_A_PCM*A_hex*k_PCM); // thermal resistance

        b_PCM= (c_PCM_m-c_s_PCM)/(T_PCM_m-T_m1); // T-dependent coefficient in c_PCM calculation
        d_PCM= (c_PCM_m - c_l_PCM)/(T_m2-T_PCM_m); // T-dependent coefficient in c_PCM calculation

        M_dx = rho_PCM*A_hex*dx_PCM;
        V_dx = A_hex*dx_PCM; // discretized volume

        if Convection == true then
           k_eff_ratio = max(0.386*(Pr/(0.861 + Pr))^(1/4)*Ra^(1/4),1);
        k_PCM_m = k_s_PCM + (k_l_PCM*k_eff_ratio2-k_s_PCM)*(T_PCM_m - T_m1)/(T_m2 -T_m1); // linearly interpolated thermal conductivity at T = T_m
         k_eff_ratio2 = max(0.386*(Pr/(0.861 + Pr))^(1/4)*Ra2^(1/4), 1);

        else
          k_eff_ratio = 1;
           k_eff_ratio2 = 1;

        k_PCM_m = k_s_PCM + (k_l_PCM-k_s_PCM)*(T_PCM_m - T_m1)/(T_m2 -T_m1);

        end if;

        // The specific heat capacity of the PCM is calculated depending on the physical condition (i.e. solid, two-phase: melting, liquid)
        if T_PCM < T_m1 then
          c_PCM = c_s_PCM;
          rho_PCM = rho_s_PCM;
          k_PCM = k_s_PCM;
          //Ra = 0;
          SOLID = true;
          TWO_PHASE1 = false;
          TWO_PHASE2 = false;
          LIQUID = false;

        else
          if T_m1<= T_PCM and T_PCM <= T_PCM_m then
          c_PCM= a + b_PCM*(T_PCM-T_m1);
          rho_PCM = rho_s_PCM + (rho_PCM_m - rho_s_PCM)*(T_PCM - T_m1)/(T_PCM_m-T_m1);
          k_PCM = k_s_PCM + (k_PCM_m -k_s_PCM)*(T_PCM - T_m1)/(T_PCM_m-T_m1);
          //k_PCM_m
          //Ra = 0;
          SOLID = false;
          TWO_PHASE1 = true;
          TWO_PHASE2 = false;
          LIQUID = false;

        //  k_eff_ratio = 1;

        //f=h_PCM/lambda;

        else if T_PCM_m < T_PCM and T_PCM< T_m2 then
          c_PCM = c - d_PCM*(T_PCM-T_PCM_m);
          rho_PCM = rho_PCM_m + (rho_l_PCM-rho_PCM_m)*(T_PCM - T_PCM_m)/(T_m2-T_PCM_m);
          k_PCM = k_PCM_m + ((k_l_PCM*k_eff_ratio2)-k_PCM_m)*(T_PCM - T_PCM_m)/(T_m2-T_PCM_m);

          //Ra = 0;
        //  k_eff_ratio = 1;

          SOLID = false;
          TWO_PHASE1 = false;
          TWO_PHASE2 = true;
          LIQUID = false;

        //f=h_PCM/lambda;

        else
          c_PCM = c_l_PCM;
          rho_PCM = rho_l_PCM;

        k_PCM = k_l_PCM*k_eff_ratio;

        //f=1;

          SOLID = false;
          TWO_PHASE1 = false;
          TWO_PHASE2 = false;
          LIQUID = true;

        end if;
        end if;
        end if;

         Ra =(g*beta*abs(delta_T2)*Lc^3)/(nu*k_l_PCM/(rho_l_PCM*c_l_PCM));
        // k_eff_ratio = 0.386*(Pr/(0.861 + Pr))^(1/4)*Ra^(1/4);

        k_Test = k_l_PCM*k_eff_ratio2;
        k_Test1 = k_l_PCM*k_eff_ratio;

         Ra11 =(g*beta*abs(Delta_T1)*Lc^3)/(nu*k_l_PCM/(rho_l_PCM*c_l_PCM));
         k_eff_ratio11 = max(0.386*(Pr/(0.861 + Pr))^(1/4)*Ra11^(1/4), 1);

        delta_T2 = Delta_T1;

         Ra22 =(g*beta*abs(Delta_T22)*Lc^3)/(nu*k_l_PCM/(rho_l_PCM*c_l_PCM));
         k_eff_ratio22 = 0.386*(Pr/(0.861 + Pr))^(1/4)*Ra22^(1/4);

         //Ra2 =(g*beta*(T_m2-T_PCM_start)*Lc^3)/(nu*k_l_PCM/(rho_l_PCM*c_l_PCM));
          Ra2 =max((g*beta*(Delta_T22)*Lc^3)/(nu*k_l_PCM/(rho_l_PCM*c_l_PCM)), 0);
        /* Ra33 =max((g*beta*(Delta_T33)*Lc^3)/(nu*k_l_PCM/(rho_l_PCM*c_l_PCM)),0); // to be compared with the correct one (i.e. Ra2)
 k_eff_ratio33 = max(0.386*(Pr/(0.861 + Pr))^(1/4)*Ra33^(1/4), 1); // to be compared with the correct one (i.e. k_eff_ratio_2)
*/
        // Counter is an Integer variable that keeps track of the actual position of the PCM block to calculate the correct heat transfer area
        // It is initialized with a value of 0 only for the "A" port of "PCM_inner"
        A.Counter + 1 = B.Counter;
        A.P = B.P;

        /*
[1] Farid. M.M. et al. "Melting and solidifcation in multi.dimensional geometry and presence of more than one interface".
    Energy and Conversion Management, Vol. 38, No. 8 pp. 809-818, 1998. DOI: 10.1016/S0196-8904(97)00038-1.
    http://findit.dtu.dk/en/catalog/13046605 

[2] Incropera F. et al. "Introduction to Heat transfer". Wiley&Sons 6th Ed., pp. 592.593. 
    ISBN 13 978-0470-50196-2.
*/

          annotation (Diagram(graphics), Icon(graphics={Bitmap(extent={{-70,-50},
                      {72,70}}, fileName=
                      "modelica://HySDeP/Graphics/PCM - rubrick cube figure for Dymola.jpg")}));
        end PCM;

        model Aluminum_thickness
          import SI = Modelica.SIunits;

        parameter SI.Length dx_Al=2e-3
            "assume 2mm aluminum thickness to confine the PCM";

          HeatPort A annotation (Placement(transformation(extent={{-110,-10},{-90,10}}),
                iconTransformation(extent={{-128,-26},{-90,10}})));
          HeatPort B annotation (Placement(transformation(extent={{90,-10},{110,10}}),
                iconTransformation(extent={{80,-20},{118,16}})));
        parameter SI.ThermalConductivity k_Al = 202 "W/(m K)"; // aluminium thermal conductivity
        parameter SI.Density rho_Al = 2700;
        parameter SI.SpecificHeatCapacity c_Al= 871;
        SI.ThermalResistance Rth "K/W";
        SI.Area A_Al;
        SI.Temperature T_Al(start=293.15, fixed=true);
         outer parameter SI.Length L;
         outer parameter SI.Length D;
        // parameter SI.Length L = 0.8;
        // parameter SI.Length D = 0.3785 "tank diameter in [m]";

        equation
        (A_Al*dx_Al)*rho_Al*c_Al*der(T_Al) = A.QDot + B.QDot;

        A_Al = Modelica.Constants.pi*(D+dx_Al)*L + 2*Modelica.Constants.pi*(D/2)^2;

          Rth = dx_Al/(A_Al*k_Al);

          A.QDot = (A.T-T_Al)/(Rth/2);
        B.QDot = (B.T-T_Al)/(Rth/2);

          annotation (Icon(graphics={Rectangle(
                  extent={{-20,78},{22,-58}},
                  lineColor={0,0,0},
                  fillColor={215,215,215},
                  fillPattern=FillPattern.Solid)}), Diagram(graphics));
        end Aluminum_thickness;

        model Adiabatic_Boundary

          import SI = Modelica.SIunits;

          HeatPort D
            annotation (Placement(transformation(extent={{-18,-104},{18,-72}})));

        equation
        D.QDot = 0;

          annotation (Icon(graphics={Rectangle(
                  extent={{-8,78},{12,-84}},
                  lineColor={0,0,0},
                  fillColor={0,0,0},
                  fillPattern=FillPattern.Solid)}));
        end Adiabatic_Boundary;

        function transition_factor
          "Get weighting factor for smooth transition (from 0 to 1)"
          extends Modelica.Icons.Function;
          import Modelica.Constants.pi;
          import Modelica.Constants.e;

          input Real start =    0.25 "start of transition interval";
          input Real stop =     0.75 "end of transition interval";
          input Real position = 0.0 "current position";
          input Integer order = 2 "Smooth up to which derivative?";

          output Real tFactor "weighting factor";

        protected
          Real[4] a_map = {-1/2, -2/pi, -3/4, -8/pi} "First parameters";
          Real[4] b_map = {1/2,   1/2,   1/2,  1/2} "Second parameters";

          // Rename variables to match with Richter2008, p.68ff
          Real phi "current phase";
          Real a "multiplier";
          Real b "addition";
          Real x_t "Start of transition";
          Real x "Current position";
          Real DELTAx "Length of transition";

          // Parameters for generalised logistic function
          Real A = 0 "Lower asymptote";
          Real K = 1 "Upper asymptote";
          Real B = 8 "Growth rate";
          Real nu= 1 "Symmetry changes";
          Real Q = nu "Zero correction";
          Real M = nu "Maximum growth for Q = nu";
          Real X = 0;
          Real END =     0;
          Real START =   0;
          Real factor =  0;

        algorithm
          assert(order>=0, "This function only supports positive values for the order of smooth derivatives.");
          assert(start<stop, "There is only support for positive differences, please provide start < stop.");

          // 0th to 2nd order
          a      := a_map[order+1];
          b      := b_map[order+1];
          x      := position;
          DELTAx := stop - start;
          x_t    := start + 0.5*DELTAx;
          phi    := (x - x_t) / DELTAx * pi;

          // higher order
          // We need to do some arbitrary scaling:
          END   :=  4.0;
          START := -2.0;
          factor := (END-START) / (stop-start);
          X := START + (position - start) * factor;

          tFactor := 1-smooth(5,noEvent(
          if position < start then
            1
          elseif position > stop then
            0
          else
            if (order == 0) then
              a                      * sin(phi)                                         + b
            elseif (order == 1) then
              a * ( 1/2 * cos(phi)   * sin(phi) + 1/2*phi)                              + b
            elseif (order == 2) then
              a * ( 1/3 * cos(phi)^2 * sin(phi) + 2/3 * sin(phi))                       + b
            else
              1 - (A + ( K-A)  / ( 1 + Q * e^(-B*(X - M)))^(1/nu))));
        //     elseif (order == 3) then
        //       a * ( 1/4 * cos(phi)^3 * sin(phi) + 3/8 * cos(phi) * sin(phi) + 3/8*phi)  + b

          annotation (smoothOrder=5,Documentation(info="<html>
<p><h4><font color=\"#008000\">DESCRIPTION:</font></h4></p>
<p>This function returns a value between 1 and 0. A smooth transition is achieved by means of defining the <code>position</code> and the transition intervall from <code>start</code> to <code>stop</code> parameter. Outside this intervall, the 1 and the 0 are returned, respectively. This transition function with up to two smooth derivatives was taken from [1]. If you provide an <code>order</code> higher than 2, the generalised logistic function[2] will be used to calculated the transition curve.</p>
<p><h4><font color=\"#008000\">OUTPUT:</font></h4></p>
<p><code>tFactor</code> = smooth transition between 0 and 1 from <code>start</code> to <code>stop</code> [-]</p>
<p>Use <code>tFactor</code> in an equation like this: </p>
<pre>tFactor&nbsp;=&nbsp;transition_factor(start=start,stop=stop,position=position);
smoothed = tFactor*1stValue&nbsp;+&nbsp;(1&nbsp;-&nbsp;tFactor)*2ndValue;</pre>
<p><h4><font color=\"#008000\">REFERENCES:</font></h4></p>
<p>[1] Christoph C Richter, Proposal of New Object-Oriented Equation-Based Model Libraries for Thermodynamic Systems, PhD thesis, Technical University Carolo-Wilhelmina Braunschweig, 2008 </p>
<p>[2] Generalised logistic function on <a href=\"http://en.wikipedia.org/wiki/Generalised_logistic_function\">Wikipedia</a></p>
<p><h4><font color=\"#008000\">IMPLEMENTATION:</font></h4></p>
<p>Implemented in 2012 for Technical University of Denmark, DTU Mechanical Engineering, Kongens Lyngby, Denmark by Jorrit Wronski (jowr@mek.dtu.dk)</p>
</html>",         revisions=""));
        end transition_factor;

        annotation ();
      end PCM_NEW_with_Nat_Conv_Coupling_Area;

    end HeatTransfer;

    package TESTS
      model Test_Heat_Transfer

        import SI = Modelica.SIunits;
        Tanks.Tank1 tank1_1(V=0.172)
          annotation (Placement(transformation(extent={{-42,32},{-4,46}})));
        HySDeP.CHG.SAE_J2601 hRSInfo
          annotation (Placement(transformation(extent={{22,64},{42,84}})));
        HySDeP.CHG.Sources.H2Source source1_1(T=293.15)
          annotation (Placement(transformation(extent={{20,14},{40,34}})));
        HeatTransfer.HeatTransferTank heatTransferTank(Adiabatic_Wall=true)
          annotation (Placement(transformation(extent={{-42,-28},{-4,0}})));
       replaceable package Medium =
          ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
          Modelica.Media.Interfaces.PartialMedium;
      Medium.ThermodynamicState medium;

      // Devi definire le pressioni (i.e. P_amb, ..start etc) e APRR. Inoltre risolvi il problema con T_amb (=HRSInfo.Tamb)
      inner parameter SI.Temperature  T_amb = HRSinfo.T_amb;
      inner SI.Pressure  P_amb;
      inner SI.SpecificEntropy s_0;
      inner SI.SpecificEnthalpy h_0;
      inner SI.SpecificInternalEnergy u_0;
      inner parameter SI.Pressure P_start=HRSinfo.P_start;
      inner SI.Pressure P_ref;
      inner Real APRR;
      inner SI.Pressure FP = HRSinfo.FP;
         HySDeP.CHG.SAE_J2601 HRSinfo(Fueling_protocol=1, P_start=2000000);

      equation
        connect(source1_1.portA, tank1_1.portA) annotation (Line(
            points={{30.4,21.8},{6.4,21.8},{6.4,39},{-19.7182,39}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_1.heatFlow, heatTransferTank.heatFlow) annotation (Line(
            points={{-24.7273,31.825},{-24.7273,-1.4},{-26.04,-1.4}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
       medium=Medium.setState_pT(P_amb, T_amb);

      s_0=medium.s;
      h_0=medium.h;
      u_0=h_0-P_amb*1/medium.d;

      HRSinfo.P_amb=P_amb;
      //HRSinfo.P_start=P_start;
      HRSinfo.APRR=APRR;
      HRSinfo.P_ref=P_ref;

      //HRSinfo.T_amb=T_amb;

        annotation (Diagram(graphics));
      end Test_Heat_Transfer;

      model Test_Heat_Transfer_Area_middle

        import SI = Modelica.SIunits;
        HySDeP.CHG.SAE_J2601 hRSInfo
          annotation (Placement(transformation(extent={{40,66},{66,96}})));
        HySDeP.CHG.Sources.H2Source source1_1(T=293.15)
          annotation (Placement(transformation(extent={{-40,56},{14,96}})));
       replaceable package Medium =
          ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
          Modelica.Media.Interfaces.PartialMedium;
      Medium.ThermodynamicState medium;

      parameter SI.Length dInner=0.4 "Inside diameter of cylinder" annotation(Dialog(group= "Geometry"));
       parameter SI.Length LInner = 1 "Inside length of tank" annotation(Dialog(group= "Geometry"));
       parameter Boolean Area = true
          "If true, area is calcualted from length and diameter" annotation(Dialog(group= "Geometry"));
       parameter SI.Area AInner= 1 "Inside tank area "
                                    annotation(Dialog(group= "Geometry", enable = Area == false));

      inner parameter SI.Temperature  T_amb = HRSinfo.T_amb;
      inner SI.Pressure  P_amb;
      inner SI.SpecificEntropy s_0;
      inner SI.SpecificEnthalpy h_0;
      inner SI.SpecificInternalEnergy u_0;
      inner parameter SI.Pressure P_start=HRSinfo.P_start;
      inner SI.Pressure P_ref;
      inner Real APRR;
      inner SI.Pressure FP = HRSinfo.FP;
         HySDeP.CHG.SAE_J2601 HRSinfo(Fueling_protocol=1, P_start=2000000);

      /*********************************************************************************************************/
      /***********************************   GEOMETRY   **********************************************/
      /*********************************************************************************************************/

      //inner SI.Length   d = (if Area == false then AInner/(Modelica.Constants.pi) else dInner);
      inner SI.Length   d = (if Area == false then AInner/(Modelica.Constants.pi*L + Modelica.Constants.pi*d/2) else dInner);

       inner SI.Length   L= (if Area == false then 1 else LInner);
       inner SI.Area   A = (if Area == false then AInner else 2*Modelica.Constants.pi*d/2*L+2*(d/2)^2*Modelica.Constants.pi);
      inner SI.Volume V = Modelica.Constants.pi*(d/2)^2*L;

        HySDeP.CHG.HeatTransfer.HeatTransferTank_no_PCM
                                                  heatTransferTank_Area_middle(
            Adiabatic_Wall=true,
          tank=4,
          h_charging=150)
          annotation (Placement(transformation(extent={{-50,-20},{-14,14}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation(Adiabatic=false)
          annotation (Placement(transformation(extent={{-60,32},{-28,44}})));
      equation
       medium=Medium.setState_pT(P_amb, T_amb);

      s_0=medium.s;
      h_0=medium.h;
      u_0=h_0-P_amb*1/medium.d;

      HRSinfo.P_amb=P_amb;
      //HRSinfo.P_start=P_start;
      HRSinfo.APRR=APRR;
      HRSinfo.P_ref=P_ref;

      //HRSinfo.T_amb=T_amb;

        connect(tank1_VolumeCalculation.heatFlow, heatTransferTank_Area_middle.heatFlow)
          annotation (Line(
            points={{-45.4545,32.15},{-45.4545,22.075},{-34.88,22.075},{-34.88,
                12.3}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));

        connect(tank1_VolumeCalculation.portA, source1_1.portA) annotation (Line(
            points={{-41.5273,38},{-12,38},{-12,71.6},{-11.92,71.6}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(extent={{-100,-100},{100,100}}),
                            graphics), Icon(coordinateSystem(extent={{-100,-100},
                  {100,100}})));
      end Test_Heat_Transfer_Area_middle;

      model Test_Heat_Transfer_PCM_Area_middle

        import SI = Modelica.SIunits;
        HySDeP.CHG.SAE_J2601 hRSInfo
          annotation (Placement(transformation(extent={{32,-22},{102,114}})));
        HySDeP.CHG.Sources.H2Source source1_1(T=293.15)
          annotation (Placement(transformation(extent={{-70,-26},{40,78}})));
       replaceable package Medium =
          ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
          Modelica.Media.Interfaces.PartialMedium;
      Medium.ThermodynamicState medium;

      parameter SI.Length dInner=0.4 "Inside diameter of cylinder" annotation(Dialog(group= "Geometry"));
       parameter SI.Length LInner = 1 "Inside length of tank" annotation(Dialog(group= "Geometry"));
       parameter Boolean Area = true
          "If true, area is calcualted from length and diameter" annotation(Dialog(group= "Geometry"));
       parameter SI.Area AInner= 1 "Inside tank area "    annotation(Dialog(group= "Geometry", enable = Area == false));

      inner parameter SI.Temperature  T_amb = HRSinfo.T_amb;
      inner SI.Pressure  P_amb;
      inner SI.SpecificEntropy s_0;
      inner SI.SpecificEnthalpy h_0;
      inner SI.SpecificInternalEnergy u_0;
      inner parameter SI.Pressure P_start=HRSinfo.P_start;
      inner SI.Pressure P_ref;
      inner Real APRR;
      inner SI.Pressure FP = HRSinfo.FP;
         HySDeP.CHG.SAE_J2601 HRSinfo(Fueling_protocol=1, P_start=2000000);
      //inner parameter SI.Temperature T = 293.15;

      inner SI.Length   d = (if Area == false then AInner/(Modelica.Constants.pi*L + Modelica.Constants.pi*d/2) else dInner);
       inner SI.Length   L= (if Area == false then 1 else LInner);
       inner SI.Area   A = (if Area == false then AInner else 2*Modelica.Constants.pi*d/2*L+2*(d/2)^2*Modelica.Constants.pi);
      inner SI.Volume V = Modelica.Constants.pi*(d/2)^2*L;

        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation(Adiabatic=false)
          annotation (Placement(transformation(extent={{-44,-36},{6,-18}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle(
          G2=1,
          Manual_Area_inner=true,
          Manual_Area_PCM=true,
          Go=1,
          G_Inner=15,
          G_PCM=15)
          annotation (Placement(transformation(extent={{-50,-114},{18,-42}})));
      equation
       medium=Medium.setState_pT(P_amb, T_amb);

      s_0=medium.s;
      h_0=medium.h;
      u_0=h_0-P_amb*1/medium.d;

      HRSinfo.P_amb=P_amb;
      //HRSinfo.P_start=P_start;
      HRSinfo.APRR=APRR;
      HRSinfo.P_ref=P_ref;

      //HRSinfo.T_amb=T_amb;

        connect(tank1_VolumeCalculation.portA, source1_1.portA) annotation (Line(
            points={{-15.1364,-27},{-11.8836,-27},{-11.8836,14.56},{-12.8,14.56}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(heatTransferTank_PCM_Area_middle.heatFlow,
          tank1_VolumeCalculation.heatFlow) annotation (Line(
            points={{-21.44,-45.6},{-21.44,-35.775},{-21.2727,-35.775}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (Diagram(graphics),
          experiment(
            StopTime=10000,
            NumberOfIntervals=10000,
            Tolerance=1e-006),
          __Dymola_experimentSetupOutput);
      end Test_Heat_Transfer_PCM_Area_middle;

      model TESTS_Conductivity

        import SI = Modelica.SIunits;
        HySDeP.CHG.SAE_J2601 hRSInfo
          annotation (Placement(transformation(extent={{18,80},{30,94}})));
        HySDeP.CHG.Sources.H2Source source1_1(T=293.15)
          annotation (Placement(transformation(extent={{-2,82},{14,100}})));
       replaceable package Medium =
          ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
          Modelica.Media.Interfaces.PartialMedium;
      Medium.ThermodynamicState medium;

      parameter SI.Length dInner=0.4 "Inside diameter of cylinder" annotation(Dialog(group= "Geometry"));
       parameter SI.Length LInner = 1 "Inside length of tank" annotation(Dialog(group= "Geometry"));
       parameter Boolean Area = true
          "If true, area is calcualted from length and diameter" annotation(Dialog(group= "Geometry"));
       parameter SI.Area AInner= 1 "Inside tank area "    annotation(Dialog(group= "Geometry", enable = Area == false));

      inner parameter SI.Temperature  T_amb = HRSinfo.T_amb;
      inner SI.Pressure  P_amb;
      inner SI.SpecificEntropy s_0;
      inner SI.SpecificEnthalpy h_0;
      inner SI.SpecificInternalEnergy u_0;
      inner parameter SI.Pressure P_start=HRSinfo.P_start;
      inner SI.Pressure P_ref;
      inner Real APRR;
      inner SI.Pressure FP = HRSinfo.FP;
         HySDeP.CHG.SAE_J2601 HRSinfo(Fueling_protocol=1, P_start=2000000);
      //inner parameter SI.Temperature T = 293.15;

      inner SI.Length   d = (if Area == false then AInner/(Modelica.Constants.pi*L + Modelica.Constants.pi*d/2) else dInner);
       inner SI.Length   L= (if Area == false then 1 else LInner);
       inner SI.Area   A = (if Area == false then AInner else 2*Modelica.Constants.pi*d/2*L+2*(d/2)^2*Modelica.Constants.pi);
      inner SI.Volume V = Modelica.Constants.pi*(d/2)^2*L;

        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle(
          h_charging=150,
          N=100,
          t_PCM=1e-2,
          Go=1)
          annotation (Placement(transformation(extent={{-68,60},{-44,86}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation(Adiabatic=false)
          annotation (Placement(transformation(extent={{-56,86},{-26,96}})));
        HySDeP.CHG.SAE_J2601 hRSInfo1
          annotation (Placement(transformation(extent={{16,40},{28,54}})));
        HySDeP.CHG.Sources.H2Source source1_2(T=293.15)
          annotation (Placement(transformation(extent={{-4,42},{12,60}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle1(
          h_charging=150,
          N=100,
          t_PCM=1e-2,
          Go=2)
          annotation (Placement(transformation(extent={{-70,20},{-46,46}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation1(Adiabatic=false)
          annotation (Placement(transformation(extent={{-58,46},{-28,56}})));
        HySDeP.CHG.SAE_J2601 hRSInfo2
          annotation (Placement(transformation(extent={{14,-4},{26,10}})));
        HySDeP.CHG.Sources.H2Source source1_3(T=293.15)
          annotation (Placement(transformation(extent={{-6,-2},{10,16}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle2(
          h_charging=150,
          N=100,
          t_PCM=1e-2,
          Go=10)
          annotation (Placement(transformation(extent={{-72,-26},{-48,0}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation2(Adiabatic=false)
          annotation (Placement(transformation(extent={{-60,2},{-30,12}})));
        HySDeP.CHG.SAE_J2601 hRSInfo3
          annotation (Placement(transformation(extent={{12,-46},{24,-32}})));
        HySDeP.CHG.Sources.H2Source source1_4(T=293.15)
          annotation (Placement(transformation(extent={{-8,-44},{8,-26}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle3(
          h_charging=150,
          N=100,
          t_PCM=1e-2,
          Go=100)
          annotation (Placement(transformation(extent={{-74,-66},{-50,-40}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation3(Adiabatic=false)
          annotation (Placement(transformation(extent={{-62,-40},{-32,-30}})));
        HySDeP.CHG.SAE_J2601 hRSInfo4
          annotation (Placement(transformation(extent={{14,-86},{26,-72}})));
        HySDeP.CHG.Sources.H2Source source1_5(T=293.15)
          annotation (Placement(transformation(extent={{-6,-84},{10,-66}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle4(
          h_charging=150,
          N=100,
          t_PCM=1e-2,
          Go=1e3)
          annotation (Placement(transformation(extent={{-72,-106},{-48,-80}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation4(Adiabatic=false)
          annotation (Placement(transformation(extent={{-60,-80},{-30,-70}})));
      equation
       medium=Medium.setState_pT(P_amb, T_amb);

      s_0=medium.s;
      h_0=medium.h;
      u_0=h_0-P_amb*1/medium.d;

      HRSinfo.P_amb=P_amb;
      //HRSinfo.P_start=P_start;
      HRSinfo.APRR=APRR;
      HRSinfo.P_ref=P_ref;

      //HRSinfo.T_amb=T_amb;

        connect(tank1_VolumeCalculation.heatFlow, heatTransferTank_PCM_Area_middle.heatFlow)
          annotation (Line(
            points={{-42.3636,86.125},{-42.3636,83.075},{-57.92,83.075},{-57.92,
                84.7}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));

        connect(tank1_VolumeCalculation.portA, source1_1.portA) annotation (Line(
            points={{-38.6818,91},{-11.8836,91},{-11.8836,91.18},{-1.36,91.18}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation1.heatFlow,
          heatTransferTank_PCM_Area_middle1.heatFlow)
          annotation (Line(
            points={{-44.3636,46.125},{-44.3636,43.075},{-59.92,43.075},{-59.92,
                44.7}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation1.portA, source1_2.portA)
                                                                annotation (Line(
            points={{-40.6818,51},{-13.8836,51},{-13.8836,51.18},{-3.36,51.18}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation2.heatFlow,
          heatTransferTank_PCM_Area_middle2.heatFlow)
          annotation (Line(
            points={{-46.3636,2.125},{-46.3636,-0.925},{-61.92,-0.925},{-61.92,
                -1.3}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation2.portA, source1_3.portA)
                                                                annotation (Line(
            points={{-42.6818,7},{-15.8836,7},{-15.8836,7.18},{-5.36,7.18}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation3.heatFlow,
          heatTransferTank_PCM_Area_middle3.heatFlow)
          annotation (Line(
            points={{-48.3636,-39.875},{-48.3636,-42.925},{-63.92,-42.925},{
                -63.92,-41.3}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation3.portA, source1_4.portA)
                                                                annotation (Line(
            points={{-44.6818,-35},{-17.8836,-35},{-17.8836,-34.82},{-7.36,
                -34.82}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation4.heatFlow,
          heatTransferTank_PCM_Area_middle4.heatFlow)
          annotation (Line(
            points={{-46.3636,-79.875},{-46.3636,-82.925},{-61.92,-82.925},{
                -61.92,-81.3}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation4.portA, source1_5.portA)
                                                                annotation (Line(
            points={{-42.6818,-75},{-15.8836,-75},{-15.8836,-74.82},{-5.36,
                -74.82}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (Diagram(graphics),
          experiment(
            StopTime=10000,
            NumberOfIntervals=10000,
            Tolerance=1e-006),
          __Dymola_experimentSetupOutput);
      end TESTS_Conductivity;

      model TESTS_Hydrogen_Convection_Coefficient

        import SI = Modelica.SIunits;
        HySDeP.CHG.SAE_J2601 hRSInfo
          annotation (Placement(transformation(extent={{18,80},{30,94}})));
        HySDeP.CHG.Sources.H2Source source1_1(T=293.15)
          annotation (Placement(transformation(extent={{-2,82},{14,100}})));
       replaceable package Medium =
          ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
          Modelica.Media.Interfaces.PartialMedium;
      Medium.ThermodynamicState medium;

      parameter SI.Length dInner=0.4 "Inside diameter of cylinder" annotation(Dialog(group= "Geometry"));
       parameter SI.Length LInner = 1 "Inside length of tank" annotation(Dialog(group= "Geometry"));
       parameter Boolean Area = true
          "If true, area is calcualted from length and diameter" annotation(Dialog(group= "Geometry"));
       parameter SI.Area AInner= 1 "Inside tank area "    annotation(Dialog(group= "Geometry", enable = Area == false));

      inner parameter SI.Temperature  T_amb = HRSinfo.T_amb;
      inner SI.Pressure  P_amb;
      inner SI.SpecificEntropy s_0;
      inner SI.SpecificEnthalpy h_0;
      inner SI.SpecificInternalEnergy u_0;
      inner parameter SI.Pressure P_start=HRSinfo.P_start;
      inner SI.Pressure P_ref;
      inner Real APRR;
      inner SI.Pressure FP = HRSinfo.FP;
         HySDeP.CHG.SAE_J2601 HRSinfo(Fueling_protocol=1, P_start=2000000);
      //inner parameter SI.Temperature T = 293.15;

      inner SI.Length   d = (if Area == false then AInner/(Modelica.Constants.pi*L + Modelica.Constants.pi*d/2) else dInner);
       inner SI.Length   L= (if Area == false then 1 else LInner);
       inner SI.Area   A = (if Area == false then AInner else 2*Modelica.Constants.pi*d/2*L+2*(d/2)^2*Modelica.Constants.pi);
      inner SI.Volume V = Modelica.Constants.pi*(d/2)^2*L;

        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle(
          h_charging=150,
          N=100,
          t_PCM=1e-2,
          Go=1)
          annotation (Placement(transformation(extent={{-68,60},{-44,86}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation(Adiabatic=false)
          annotation (Placement(transformation(extent={{-56,86},{-26,96}})));
        HySDeP.CHG.SAE_J2601 hRSInfo1
          annotation (Placement(transformation(extent={{16,40},{28,54}})));
        HySDeP.CHG.Sources.H2Source source1_2(T=293.15)
          annotation (Placement(transformation(extent={{-4,42},{12,60}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle1(
          N=100,
          t_PCM=1e-2,
          h_charging=500,
          Go=1)
          annotation (Placement(transformation(extent={{-70,18},{-46,44}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation1(Adiabatic=false)
          annotation (Placement(transformation(extent={{-58,46},{-28,56}})));
        HySDeP.CHG.SAE_J2601 hRSInfo2
          annotation (Placement(transformation(extent={{14,-4},{26,10}})));
        HySDeP.CHG.Sources.H2Source source1_3(T=293.15)
          annotation (Placement(transformation(extent={{-6,-2},{10,16}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle2(
          N=100,
          t_PCM=1e-2,
          h_charging=1500,
          Go=1)
          annotation (Placement(transformation(extent={{-72,-24},{-48,2}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation2(Adiabatic=false)
          annotation (Placement(transformation(extent={{-60,2},{-30,12}})));
      equation
       medium=Medium.setState_pT(P_amb, T_amb);

      s_0=medium.s;
      h_0=medium.h;
      u_0=h_0-P_amb*1/medium.d;

      HRSinfo.P_amb=P_amb;
      //HRSinfo.P_start=P_start;
      HRSinfo.APRR=APRR;
      HRSinfo.P_ref=P_ref;

      //HRSinfo.T_amb=T_amb;

        connect(tank1_VolumeCalculation.heatFlow, heatTransferTank_PCM_Area_middle.heatFlow)
          annotation (Line(
            points={{-42.3636,86.125},{-42.3636,83.075},{-57.92,83.075},{-57.92,
                84.7}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));

        connect(tank1_VolumeCalculation.portA, source1_1.portA) annotation (Line(
            points={{-38.6818,91},{-11.8836,91},{-11.8836,89.02},{6.32,89.02}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation1.heatFlow,
          heatTransferTank_PCM_Area_middle1.heatFlow)
          annotation (Line(
            points={{-44.3636,46.125},{-44.3636,43.075},{-59.92,43.075},{-59.92,
                42.7}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation1.portA, source1_2.portA)
                                                                annotation (Line(
            points={{-40.6818,51},{-13.8836,51},{-13.8836,49.02},{4.32,49.02}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation2.heatFlow,
          heatTransferTank_PCM_Area_middle2.heatFlow)
          annotation (Line(
            points={{-46.3636,2.125},{-46.3636,-0.925},{-61.92,-0.925},{-61.92,
                0.7}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation2.portA, source1_3.portA)
                                                                annotation (Line(
            points={{-42.6818,7},{-15.8836,7},{-15.8836,5.02},{2.32,5.02}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (Diagram(graphics),
          experiment(
            StopTime=10000,
            NumberOfIntervals=10000,
            Tolerance=1e-006),
          __Dymola_experimentSetupOutput);
      end TESTS_Hydrogen_Convection_Coefficient;

      model TESTS_Thickness

        import SI = Modelica.SIunits;
        HySDeP.CHG.SAE_J2601 hRSInfo
          annotation (Placement(transformation(extent={{18,80},{30,94}})));
        HySDeP.CHG.Sources.H2Source source1_1(T=293.15)
          annotation (Placement(transformation(extent={{-2,82},{14,100}})));
       replaceable package Medium =
          ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
          Modelica.Media.Interfaces.PartialMedium;
      Medium.ThermodynamicState medium;

      parameter SI.Length dInner=0.4 "Inside diameter of cylinder" annotation(Dialog(group= "Geometry"));
       parameter SI.Length LInner = 1 "Inside length of tank" annotation(Dialog(group= "Geometry"));
       parameter Boolean Area = true
          "If true, area is calcualted from length and diameter" annotation(Dialog(group= "Geometry"));
       parameter SI.Area AInner= 1 "Inside tank area "    annotation(Dialog(group= "Geometry", enable = Area == false));

      inner parameter SI.Temperature  T_amb = HRSinfo.T_amb;
      inner SI.Pressure  P_amb;
      inner SI.SpecificEntropy s_0;
      inner SI.SpecificEnthalpy h_0;
      inner SI.SpecificInternalEnergy u_0;
      inner parameter SI.Pressure P_start=HRSinfo.P_start;
      inner SI.Pressure P_ref;
      inner Real APRR;
      inner SI.Pressure FP = HRSinfo.FP;
         HySDeP.CHG.SAE_J2601 HRSinfo(Fueling_protocol=1, P_start=2000000);
      //inner parameter SI.Temperature T = 293.15;

      inner SI.Length   d = (if Area == false then AInner/(Modelica.Constants.pi*L + Modelica.Constants.pi*d/2) else dInner);
       inner SI.Length   L= (if Area == false then 1 else LInner);
       inner SI.Area   A = (if Area == false then AInner else 2*Modelica.Constants.pi*d/2*L+2*(d/2)^2*Modelica.Constants.pi);
      inner SI.Volume V = Modelica.Constants.pi*(d/2)^2*L;

        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle(
          h_charging=150,
          N=100,
          t_PCM=1e-2,
          Go=1)
          annotation (Placement(transformation(extent={{-68,60},{-44,86}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation(Adiabatic=false)
          annotation (Placement(transformation(extent={{-56,86},{-26,96}})));
        HySDeP.CHG.SAE_J2601 hRSInfo1
          annotation (Placement(transformation(extent={{16,40},{28,54}})));
        HySDeP.CHG.Sources.H2Source source1_2(T=293.15)
          annotation (Placement(transformation(extent={{-4,42},{12,60}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle1(
          h_charging=150,
          Go=1,
          N=200,
          t_PCM=2e-2)
          annotation (Placement(transformation(extent={{-70,20},{-46,46}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation1(Adiabatic=false)
          annotation (Placement(transformation(extent={{-58,46},{-28,56}})));
        HySDeP.CHG.SAE_J2601 hRSInfo2
          annotation (Placement(transformation(extent={{14,-4},{26,10}})));
        HySDeP.CHG.Sources.H2Source source1_3(T=293.15)
          annotation (Placement(transformation(extent={{-6,-2},{10,16}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle2(
          h_charging=150,
          Go=1,
          N=500,
          t_PCM=5e-2)
          annotation (Placement(transformation(extent={{-72,-24},{-48,2}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation2(Adiabatic=false)
          annotation (Placement(transformation(extent={{-60,2},{-30,12}})));
        HySDeP.CHG.SAE_J2601 hRSInfo3
          annotation (Placement(transformation(extent={{12,-46},{24,-32}})));
        HySDeP.CHG.Sources.H2Source source1_4(T=293.15)
          annotation (Placement(transformation(extent={{-8,-44},{8,-26}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle3(
          h_charging=150,
          Go=1,
          N=500,
          t_PCM=10e-2)
          annotation (Placement(transformation(extent={{-74,-66},{-50,-40}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation3(Adiabatic=false)
          annotation (Placement(transformation(extent={{-62,-40},{-32,-30}})));
        HySDeP.CHG.SAE_J2601 hRSInfo4
          annotation (Placement(transformation(extent={{14,-86},{26,-72}})));
        HySDeP.CHG.Sources.H2Source source1_5(T=293.15)
          annotation (Placement(transformation(extent={{-6,-84},{10,-66}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle4(
          h_charging=150,
          Go=1,
          N=500,
          t_PCM=20e-2)
          annotation (Placement(transformation(extent={{-72,-106},{-48,-80}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation4(Adiabatic=false)
          annotation (Placement(transformation(extent={{-60,-80},{-30,-70}})));
      equation
       medium=Medium.setState_pT(P_amb, T_amb);

      s_0=medium.s;
      h_0=medium.h;
      u_0=h_0-P_amb*1/medium.d;

      HRSinfo.P_amb=P_amb;
      //HRSinfo.P_start=P_start;
      HRSinfo.APRR=APRR;
      HRSinfo.P_ref=P_ref;

      //HRSinfo.T_amb=T_amb;

        connect(tank1_VolumeCalculation.heatFlow, heatTransferTank_PCM_Area_middle.heatFlow)
          annotation (Line(
            points={{-42.3636,86.125},{-42.3636,83.075},{-57.92,83.075},{-57.92,
                84.7}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));

        connect(tank1_VolumeCalculation.portA, source1_1.portA) annotation (Line(
            points={{-38.6818,91},{-11.8836,91},{-11.8836,91.18},{-1.36,91.18}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation1.heatFlow,
          heatTransferTank_PCM_Area_middle1.heatFlow)
          annotation (Line(
            points={{-44.3636,46.125},{-44.3636,43.075},{-59.92,43.075},{-59.92,
                44.7}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation1.portA, source1_2.portA)
                                                                annotation (Line(
            points={{-40.6818,51},{-13.8836,51},{-13.8836,51.18},{-3.36,51.18}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation2.heatFlow,
          heatTransferTank_PCM_Area_middle2.heatFlow)
          annotation (Line(
            points={{-46.3636,2.125},{-46.3636,-0.925},{-61.92,-0.925},{-61.92,
                0.7}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation2.portA, source1_3.portA)
                                                                annotation (Line(
            points={{-42.6818,7},{-15.8836,7},{-15.8836,7.18},{-5.36,7.18}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation3.heatFlow,
          heatTransferTank_PCM_Area_middle3.heatFlow)
          annotation (Line(
            points={{-48.3636,-39.875},{-48.3636,-42.925},{-63.92,-42.925},{
                -63.92,-41.3}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation3.portA, source1_4.portA)
                                                                annotation (Line(
            points={{-44.6818,-35},{-17.8836,-35},{-17.8836,-34.82},{-7.36,
                -34.82}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation4.heatFlow,
          heatTransferTank_PCM_Area_middle4.heatFlow)
          annotation (Line(
            points={{-46.3636,-79.875},{-46.3636,-82.925},{-61.92,-82.925},{
                -61.92,-81.3}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation4.portA, source1_5.portA)
                                                                annotation (Line(
            points={{-42.6818,-75},{-15.8836,-75},{-15.8836,-74.82},{-5.36,
                -74.82}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (Diagram(graphics),
          experiment(
            StopTime=10000,
            NumberOfIntervals=10000,
            Tolerance=1e-006),
          __Dymola_experimentSetupOutput);
      end TESTS_Thickness;

      model TESTS_density

        import SI = Modelica.SIunits;
        HySDeP.CHG.SAE_J2601 hRSInfo
          annotation (Placement(transformation(extent={{18,80},{30,94}})));
        HySDeP.CHG.Sources.H2Source source1_1(T=293.15)
          annotation (Placement(transformation(extent={{-2,82},{14,100}})));
       replaceable package Medium =
          ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
          Modelica.Media.Interfaces.PartialMedium;
      Medium.ThermodynamicState medium;

      parameter SI.Length dInner=0.4 "Inside diameter of cylinder" annotation(Dialog(group= "Geometry"));
       parameter SI.Length LInner = 1 "Inside length of tank" annotation(Dialog(group= "Geometry"));
       parameter Boolean Area = true
          "If true, area is calcualted from length and diameter" annotation(Dialog(group= "Geometry"));
       parameter SI.Area AInner= 1 "Inside tank area "    annotation(Dialog(group= "Geometry", enable = Area == false));

      inner parameter SI.Temperature  T_amb = HRSinfo.T_amb;
      inner SI.Pressure  P_amb;
      inner SI.SpecificEntropy s_0;
      inner SI.SpecificEnthalpy h_0;
      inner SI.SpecificInternalEnergy u_0;
      inner parameter SI.Pressure P_start=HRSinfo.P_start;
      inner SI.Pressure P_ref;
      inner Real APRR;
      inner SI.Pressure FP = HRSinfo.FP;
         HySDeP.CHG.SAE_J2601 HRSinfo(Fueling_protocol=1, P_start=2000000);
      //inner parameter SI.Temperature T = 293.15;

      inner SI.Length   d = (if Area == false then AInner/(Modelica.Constants.pi*L + Modelica.Constants.pi*d/2) else dInner);
       inner SI.Length   L= (if Area == false then 1 else LInner);
       inner SI.Area   A = (if Area == false then AInner else 2*Modelica.Constants.pi*d/2*L+2*(d/2)^2*Modelica.Constants.pi);
      inner SI.Volume V = Modelica.Constants.pi*(d/2)^2*L;

        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle(
          h_charging=150,
          N=100,
          t_PCM=1e-2,
          Go=1)
          annotation (Placement(transformation(extent={{-68,60},{-44,86}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation(Adiabatic=false)
          annotation (Placement(transformation(extent={{-56,86},{-26,96}})));
        HySDeP.CHG.SAE_J2601 hRSInfo1
          annotation (Placement(transformation(extent={{16,40},{28,54}})));
        HySDeP.CHG.Sources.H2Source source1_2(T=293.15)
          annotation (Placement(transformation(extent={{-4,42},{12,60}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle1(
          h_charging=150,
          N=100,
          t_PCM=1e-2,
          Go=1,
          G1=2)
          annotation (Placement(transformation(extent={{-70,20},{-46,46}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation1(Adiabatic=false)
          annotation (Placement(transformation(extent={{-58,46},{-28,56}})));
        HySDeP.CHG.SAE_J2601 hRSInfo2
          annotation (Placement(transformation(extent={{14,-4},{26,10}})));
        HySDeP.CHG.Sources.H2Source source1_3(T=293.15)
          annotation (Placement(transformation(extent={{-6,-2},{10,16}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle2(
          h_charging=150,
          N=100,
          t_PCM=1e-2,
          Go=1,
          G1=10)
          annotation (Placement(transformation(extent={{-72,-24},{-48,2}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation2(Adiabatic=false)
          annotation (Placement(transformation(extent={{-60,2},{-30,12}})));
        HySDeP.CHG.SAE_J2601 hRSInfo3
          annotation (Placement(transformation(extent={{12,-46},{24,-32}})));
        HySDeP.CHG.Sources.H2Source source1_4(T=293.15)
          annotation (Placement(transformation(extent={{-8,-44},{8,-26}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle3(
          h_charging=150,
          N=100,
          t_PCM=1e-2,
          Go=1,
          G1=100)
          annotation (Placement(transformation(extent={{-74,-66},{-50,-40}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation3(Adiabatic=false)
          annotation (Placement(transformation(extent={{-62,-40},{-32,-30}})));
      equation
       medium=Medium.setState_pT(P_amb, T_amb);

      s_0=medium.s;
      h_0=medium.h;
      u_0=h_0-P_amb*1/medium.d;

      HRSinfo.P_amb=P_amb;
      //HRSinfo.P_start=P_start;
      HRSinfo.APRR=APRR;
      HRSinfo.P_ref=P_ref;

      //HRSinfo.T_amb=T_amb;

        connect(tank1_VolumeCalculation.heatFlow, heatTransferTank_PCM_Area_middle.heatFlow)
          annotation (Line(
            points={{-42.3636,86.125},{-42.3636,83.075},{-57.92,83.075},{-57.92,
                84.7}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));

        connect(tank1_VolumeCalculation.portA, source1_1.portA) annotation (Line(
            points={{-38.6818,91},{-11.8836,91},{-11.8836,91.18},{-1.36,91.18}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation1.heatFlow,
          heatTransferTank_PCM_Area_middle1.heatFlow)
          annotation (Line(
            points={{-44.3636,46.125},{-44.3636,43.075},{-59.92,43.075},{-59.92,
                44.7}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation1.portA, source1_2.portA)
                                                                annotation (Line(
            points={{-40.6818,51},{-13.8836,51},{-13.8836,51.18},{-3.36,51.18}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation2.heatFlow,
          heatTransferTank_PCM_Area_middle2.heatFlow)
          annotation (Line(
            points={{-46.3636,2.125},{-46.3636,-0.925},{-61.92,-0.925},{-61.92,
                0.7}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation2.portA, source1_3.portA)
                                                                annotation (Line(
            points={{-42.6818,7},{-15.8836,7},{-15.8836,7.18},{-5.36,7.18}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation3.heatFlow,
          heatTransferTank_PCM_Area_middle3.heatFlow)
          annotation (Line(
            points={{-48.3636,-39.875},{-48.3636,-42.925},{-63.92,-42.925},{
                -63.92,-41.3}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation3.portA, source1_4.portA)
                                                                annotation (Line(
            points={{-44.6818,-35},{-17.8836,-35},{-17.8836,-34.82},{-7.36,
                -34.82}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (Diagram(graphics),
          experiment(
            StopTime=10000,
            NumberOfIntervals=10000,
            Tolerance=1e-006),
          __Dymola_experimentSetupOutput);
      end TESTS_density;

      model TESTS_heat_capacity

        import SI = Modelica.SIunits;
        HySDeP.CHG.SAE_J2601 hRSInfo
          annotation (Placement(transformation(extent={{18,80},{30,94}})));
        HySDeP.CHG.Sources.H2Source source1_1(T=293.15)
          annotation (Placement(transformation(extent={{-2,82},{14,100}})));
       replaceable package Medium =
          ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
          Modelica.Media.Interfaces.PartialMedium;
      Medium.ThermodynamicState medium;

      parameter SI.Length dInner=0.4 "Inside diameter of cylinder" annotation(Dialog(group= "Geometry"));
       parameter SI.Length LInner = 1 "Inside length of tank" annotation(Dialog(group= "Geometry"));
       parameter Boolean Area = true
          "If true, area is calcualted from length and diameter" annotation(Dialog(group= "Geometry"));
       parameter SI.Area AInner= 1 "Inside tank area "    annotation(Dialog(group= "Geometry", enable = Area == false));

      inner parameter SI.Temperature  T_amb = HRSinfo.T_amb;
      inner SI.Pressure  P_amb;
      inner SI.SpecificEntropy s_0;
      inner SI.SpecificEnthalpy h_0;
      inner SI.SpecificInternalEnergy u_0;
      inner parameter SI.Pressure P_start=HRSinfo.P_start;
      inner SI.Pressure P_ref;
      inner Real APRR;
      inner SI.Pressure FP = HRSinfo.FP;
         HySDeP.CHG.SAE_J2601 HRSinfo(Fueling_protocol=1, P_start=2000000);
      //inner parameter SI.Temperature T = 293.15;

      inner SI.Length   d = (if Area == false then AInner/(Modelica.Constants.pi*L + Modelica.Constants.pi*d/2) else dInner);
       inner SI.Length   L= (if Area == false then 1 else LInner);
       inner SI.Area   A = (if Area == false then AInner else 2*Modelica.Constants.pi*d/2*L+2*(d/2)^2*Modelica.Constants.pi);
      inner SI.Volume V = Modelica.Constants.pi*(d/2)^2*L;

        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle(
          h_charging=150,
          N=100,
          t_PCM=1e-2,
          Go=1,
          G2=1)
          annotation (Placement(transformation(extent={{-68,60},{-44,86}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation(Adiabatic=false)
          annotation (Placement(transformation(extent={{-56,86},{-26,96}})));
        HySDeP.CHG.SAE_J2601 hRSInfo1
          annotation (Placement(transformation(extent={{16,40},{28,54}})));
        HySDeP.CHG.Sources.H2Source source1_2(T=293.15)
          annotation (Placement(transformation(extent={{-4,42},{12,60}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle1(
          h_charging=150,
          N=100,
          t_PCM=1e-2,
          Go=1,
          G1=1,
          G2=2)
          annotation (Placement(transformation(extent={{-70,20},{-46,46}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation1(Adiabatic=false)
          annotation (Placement(transformation(extent={{-58,46},{-28,56}})));
        HySDeP.CHG.SAE_J2601 hRSInfo2
          annotation (Placement(transformation(extent={{14,-4},{26,10}})));
        HySDeP.CHG.Sources.H2Source source1_3(T=293.15)
          annotation (Placement(transformation(extent={{-6,-2},{10,16}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle2(
          h_charging=150,
          N=100,
          t_PCM=1e-2,
          Go=1,
          G1=1,
          G2=10)
          annotation (Placement(transformation(extent={{-72,-24},{-48,2}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation2(Adiabatic=false)
          annotation (Placement(transformation(extent={{-60,2},{-30,12}})));
        HySDeP.CHG.SAE_J2601 hRSInfo3
          annotation (Placement(transformation(extent={{12,-46},{24,-32}})));
        HySDeP.CHG.Sources.H2Source source1_4(T=293.15)
          annotation (Placement(transformation(extent={{-8,-44},{8,-26}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle3(
          h_charging=150,
          N=100,
          t_PCM=1e-2,
          Go=1,
          G1=1,
          G2=100)
          annotation (Placement(transformation(extent={{-74,-66},{-50,-40}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation3(Adiabatic=false)
          annotation (Placement(transformation(extent={{-62,-40},{-32,-30}})));
      equation
       medium=Medium.setState_pT(P_amb, T_amb);

      s_0=medium.s;
      h_0=medium.h;
      u_0=h_0-P_amb*1/medium.d;

      HRSinfo.P_amb=P_amb;
      //HRSinfo.P_start=P_start;
      HRSinfo.APRR=APRR;
      HRSinfo.P_ref=P_ref;

      //HRSinfo.T_amb=T_amb;

        connect(tank1_VolumeCalculation.heatFlow, heatTransferTank_PCM_Area_middle.heatFlow)
          annotation (Line(
            points={{-42.3636,86.125},{-42.3636,83.075},{-57.92,83.075},{-57.92,
                84.7}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));

        connect(tank1_VolumeCalculation.portA, source1_1.portA) annotation (Line(
            points={{-38.6818,91},{-11.8836,91},{-11.8836,91.18},{-1.36,91.18}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation1.heatFlow,
          heatTransferTank_PCM_Area_middle1.heatFlow)
          annotation (Line(
            points={{-44.3636,46.125},{-44.3636,43.075},{-59.92,43.075},{-59.92,
                44.7}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation1.portA, source1_2.portA)
                                                                annotation (Line(
            points={{-40.6818,51},{-13.8836,51},{-13.8836,51.18},{-3.36,51.18}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation2.heatFlow,
          heatTransferTank_PCM_Area_middle2.heatFlow)
          annotation (Line(
            points={{-46.3636,2.125},{-46.3636,-0.925},{-61.92,-0.925},{-61.92,
                0.7}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation2.portA, source1_3.portA)
                                                                annotation (Line(
            points={{-42.6818,7},{-15.8836,7},{-15.8836,7.18},{-5.36,7.18}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation3.heatFlow,
          heatTransferTank_PCM_Area_middle3.heatFlow)
          annotation (Line(
            points={{-48.3636,-39.875},{-48.3636,-42.925},{-63.92,-42.925},{
                -63.92,-41.3}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation3.portA, source1_4.portA)
                                                                annotation (Line(
            points={{-44.6818,-35},{-17.8836,-35},{-17.8836,-34.82},{-7.36,
                -34.82}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (Diagram(graphics),
          experiment(
            StopTime=10000,
            NumberOfIntervals=10000,
            Tolerance=1e-006),
          __Dymola_experimentSetupOutput);
      end TESTS_heat_capacity;

      model TESTS_lambda

        import SI = Modelica.SIunits;
        HySDeP.CHG.SAE_J2601 hRSInfo
          annotation (Placement(transformation(extent={{18,80},{30,94}})));
        HySDeP.CHG.Sources.H2Source source1_1(T=293.15)
          annotation (Placement(transformation(extent={{-2,82},{14,100}})));
       replaceable package Medium =
          ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
          Modelica.Media.Interfaces.PartialMedium;
      Medium.ThermodynamicState medium;

      parameter SI.Length dInner=0.4 "Inside diameter of cylinder" annotation(Dialog(group= "Geometry"));
       parameter SI.Length LInner = 1 "Inside length of tank" annotation(Dialog(group= "Geometry"));
       parameter Boolean Area = true
          "If true, area is calcualted from length and diameter" annotation(Dialog(group= "Geometry"));
       parameter SI.Area AInner= 1 "Inside tank area "    annotation(Dialog(group= "Geometry", enable = Area == false));

      inner parameter SI.Temperature  T_amb = HRSinfo.T_amb;
      inner SI.Pressure  P_amb;
      inner SI.SpecificEntropy s_0;
      inner SI.SpecificEnthalpy h_0;
      inner SI.SpecificInternalEnergy u_0;
      inner parameter SI.Pressure P_start=HRSinfo.P_start;
      inner SI.Pressure P_ref;
      inner Real APRR;
      inner SI.Pressure FP = HRSinfo.FP;
         HySDeP.CHG.SAE_J2601 HRSinfo(Fueling_protocol=1, P_start=2000000);
      //inner parameter SI.Temperature T = 293.15;

      inner SI.Length   d = (if Area == false then AInner/(Modelica.Constants.pi*L + Modelica.Constants.pi*d/2) else dInner);
       inner SI.Length   L= (if Area == false then 1 else LInner);
       inner SI.Area   A = (if Area == false then AInner else 2*Modelica.Constants.pi*d/2*L+2*(d/2)^2*Modelica.Constants.pi);
      inner SI.Volume V = Modelica.Constants.pi*(d/2)^2*L;

        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle(
          h_charging=150,
          N=100,
          t_PCM=1e-2,
          Go=1,
          G2=1)
          annotation (Placement(transformation(extent={{-68,60},{-44,86}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation(Adiabatic=false)
          annotation (Placement(transformation(extent={{-56,86},{-26,96}})));
        HySDeP.CHG.SAE_J2601 hRSInfo1
          annotation (Placement(transformation(extent={{16,40},{28,54}})));
        HySDeP.CHG.Sources.H2Source source1_2(T=293.15)
          annotation (Placement(transformation(extent={{-4,42},{12,60}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle1(
          h_charging=150,
          N=100,
          t_PCM=1e-2,
          Go=1,
          G1=1,
          G2=1,
          G3=2)
          annotation (Placement(transformation(extent={{-70,20},{-46,46}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation1(Adiabatic=false)
          annotation (Placement(transformation(extent={{-58,46},{-28,56}})));
        HySDeP.CHG.SAE_J2601 hRSInfo2
          annotation (Placement(transformation(extent={{14,-4},{26,10}})));
        HySDeP.CHG.Sources.H2Source source1_3(T=293.15)
          annotation (Placement(transformation(extent={{-6,-2},{10,16}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle2(
          h_charging=150,
          N=100,
          t_PCM=1e-2,
          Go=1,
          G1=1,
          G2=1,
          G3=10)
          annotation (Placement(transformation(extent={{-72,-24},{-48,2}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation2(Adiabatic=false)
          annotation (Placement(transformation(extent={{-60,2},{-30,12}})));
        HySDeP.CHG.SAE_J2601 hRSInfo3
          annotation (Placement(transformation(extent={{12,-46},{24,-32}})));
        HySDeP.CHG.Sources.H2Source source1_4(T=293.15)
          annotation (Placement(transformation(extent={{-8,-44},{8,-26}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle3(
          h_charging=150,
          N=100,
          t_PCM=1e-2,
          Go=1,
          G1=1,
          G2=1,
          G3=100)
          annotation (Placement(transformation(extent={{-74,-66},{-50,-40}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation3(Adiabatic=false)
          annotation (Placement(transformation(extent={{-62,-40},{-32,-30}})));
      equation
       medium=Medium.setState_pT(P_amb, T_amb);

      s_0=medium.s;
      h_0=medium.h;
      u_0=h_0-P_amb*1/medium.d;

      HRSinfo.P_amb=P_amb;
      //HRSinfo.P_start=P_start;
      HRSinfo.APRR=APRR;
      HRSinfo.P_ref=P_ref;

      //HRSinfo.T_amb=T_amb;

        connect(tank1_VolumeCalculation.heatFlow, heatTransferTank_PCM_Area_middle.heatFlow)
          annotation (Line(
            points={{-42.3636,86.125},{-42.3636,83.075},{-57.92,83.075},{-57.92,
                84.7}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));

        connect(tank1_VolumeCalculation.portA, source1_1.portA) annotation (Line(
            points={{-38.6818,91},{-11.8836,91},{-11.8836,91.18},{-1.36,91.18}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation1.heatFlow,
          heatTransferTank_PCM_Area_middle1.heatFlow)
          annotation (Line(
            points={{-44.3636,46.125},{-44.3636,43.075},{-59.92,43.075},{-59.92,
                44.7}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation1.portA, source1_2.portA)
                                                                annotation (Line(
            points={{-40.6818,51},{-13.8836,51},{-13.8836,51.18},{-3.36,51.18}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation2.heatFlow,
          heatTransferTank_PCM_Area_middle2.heatFlow)
          annotation (Line(
            points={{-46.3636,2.125},{-46.3636,-0.925},{-61.92,-0.925},{-61.92,
                0.7}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation2.portA, source1_3.portA)
                                                                annotation (Line(
            points={{-42.6818,7},{-15.8836,7},{-15.8836,7.18},{-5.36,7.18}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation3.heatFlow,
          heatTransferTank_PCM_Area_middle3.heatFlow)
          annotation (Line(
            points={{-48.3636,-39.875},{-48.3636,-42.925},{-63.92,-42.925},{
                -63.92,-41.3}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation3.portA, source1_4.portA)
                                                                annotation (Line(
            points={{-44.6818,-35},{-17.8836,-35},{-17.8836,-34.82},{-7.36,
                -34.82}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (Diagram(graphics),
          experiment(
            StopTime=10000,
            NumberOfIntervals=10000,
            Tolerance=1e-006),
          __Dymola_experimentSetupOutput);
      end TESTS_lambda;

      model TESTS_AREA_Inner

        import SI = Modelica.SIunits;
        HySDeP.CHG.SAE_J2601 hRSInfo
          annotation (Placement(transformation(extent={{14,86},{30,100}})));
        HySDeP.CHG.Sources.H2Source source1_1(T=293.15)
          annotation (Placement(transformation(extent={{-2,66},{24,94}})));
       replaceable package Medium =
          ExternalMedia.Media.CoolPropMedium(substanceNames={"Hydrogen|calc_transport=1|enable_BICUBIC=1"}) constrainedby
          Modelica.Media.Interfaces.PartialMedium;
      Medium.ThermodynamicState medium;

      parameter SI.Length dInner=0.4 "Inside diameter of cylinder" annotation(Dialog(group= "Geometry"));
       parameter SI.Length LInner = 1 "Inside length of tank" annotation(Dialog(group= "Geometry"));
       parameter Boolean Area = true
          "If true, area is calcualted from length and diameter" annotation(Dialog(group= "Geometry"));
       parameter SI.Area AInner= 1 "Inside tank area "    annotation(Dialog(group= "Geometry", enable = Area == false));

      inner parameter SI.Temperature  T_amb = HRSinfo.T_amb;
      inner SI.Pressure  P_amb;
      inner SI.SpecificEntropy s_0;
      inner SI.SpecificEnthalpy h_0;
      inner SI.SpecificInternalEnergy u_0;
      inner parameter SI.Pressure P_start=HRSinfo.P_start;
      inner SI.Pressure P_ref;
      inner Real APRR;
      inner SI.Pressure FP = HRSinfo.FP;
         HySDeP.CHG.SAE_J2601 HRSinfo(Fueling_protocol=1, P_start=2000000);
      //inner parameter SI.Temperature T = 293.15;

      inner SI.Length   d = (if Area == false then AInner/(Modelica.Constants.pi*L + Modelica.Constants.pi*d/2) else dInner);
       inner SI.Length   L= (if Area == false then 1 else LInner);
       inner SI.Area   A = (if Area == false then AInner else 2*Modelica.Constants.pi*d/2*L+2*(d/2)^2*Modelica.Constants.pi);
      inner SI.Volume V = Modelica.Constants.pi*(d/2)^2*L;

        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation(Adiabatic=false)
          annotation (Placement(transformation(extent={{-54,74},{-16,88}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle(Manual_Area_inner=true, G_Inner=1,
          Manual_Area_PCM=true,
          G_PCM=1)
          annotation (Placement(transformation(extent={{-68,52},{-46,70}})));
        HySDeP.CHG.SAE_J2601 hRSInfo1
          annotation (Placement(transformation(extent={{14,56},{30,70}})));
        HySDeP.CHG.Sources.H2Source source1_2(T=293.15)
          annotation (Placement(transformation(extent={{-2,30},{24,58}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation1(Adiabatic=false)
          annotation (Placement(transformation(extent={{-54,38},{-16,52}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle1(Manual_Area_inner=true, G_Inner=2,
          Manual_Area_PCM=true,
          G1=1,
          G2=2)
          annotation (Placement(transformation(extent={{-68,16},{-46,34}})));
        HySDeP.CHG.SAE_J2601 hRSInfo2
          annotation (Placement(transformation(extent={{16,20},{32,34}})));
        HySDeP.CHG.Sources.H2Source source1_5(T=293.15)
          annotation (Placement(transformation(extent={{2,-82},{28,-54}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation4(Adiabatic=false)
          annotation (Placement(transformation(extent={{-52,-74},{-14,-60}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle4(Manual_Area_inner=true, G_Inner=100,
          Manual_Area_PCM=true,
          G_PCM=100,
          G2=100)
          annotation (Placement(transformation(extent={{-64,-96},{-42,-78}})));
        HySDeP.CHG.Sources.H2Source source1_3(T=293.15)
          annotation (Placement(transformation(extent={{-2,-12},{24,16}})));
        HySDeP.CHG.Tanks.Tank tank1_VolumeCalculation2(Adiabatic=false)
          annotation (Placement(transformation(extent={{-54,-4},{-16,10}})));
        HySDeP.CHG.HeatTransfer.HeatTransferTank_PCM
          heatTransferTank_PCM_Area_middle2(Manual_Area_inner=true, G_Inner=10,
          Manual_Area_PCM=true,
          G_PCM=10,
          G2=10)
          annotation (Placement(transformation(extent={{-68,-26},{-46,-8}})));
        HySDeP.CHG.SAE_J2601 hRSInfo3
          annotation (Placement(transformation(extent={{16,-22},{32,-8}})));
      equation
       medium=Medium.setState_pT(P_amb, T_amb);

      s_0=medium.s;
      h_0=medium.h;
      u_0=h_0-P_amb*1/medium.d;

      HRSinfo.P_amb=P_amb;
      //HRSinfo.P_start=P_start;
      HRSinfo.APRR=APRR;
      HRSinfo.P_ref=P_ref;

      //HRSinfo.T_amb=T_amb;

        connect(tank1_VolumeCalculation.portA, source1_1.portA) annotation (Line(
            points={{-32.0636,81},{-23.8836,81},{-23.8836,76.92},{11.52,76.92}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(heatTransferTank_PCM_Area_middle.heatFlow,
          tank1_VolumeCalculation.heatFlow) annotation (Line(
            points={{-58.76,69.1},{-58.76,74.175},{-36.7273,74.175}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation1.portA, source1_2.portA)
                                                                annotation (Line(
            points={{-32.0636,45},{-23.8836,45},{-23.8836,40.92},{11.52,40.92}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(heatTransferTank_PCM_Area_middle1.heatFlow,
          tank1_VolumeCalculation1.heatFlow) annotation (Line(
            points={{-58.76,33.1},{-58.76,38.175},{-36.7273,38.175}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation4.portA, source1_5.portA)
                                                                annotation (Line(
            points={{-30.0636,-67},{-19.8836,-67},{-19.8836,-71.08},{15.52,
                -71.08}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(heatTransferTank_PCM_Area_middle4.heatFlow,
          tank1_VolumeCalculation4.heatFlow) annotation (Line(
            points={{-54.76,-78.9},{-54.76,-73.825},{-34.7273,-73.825}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(tank1_VolumeCalculation2.portA,source1_3. portA)
                                                                annotation (Line(
            points={{-32.0636,3},{-23.8836,3},{-23.8836,-1.08},{11.52,-1.08}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        connect(heatTransferTank_PCM_Area_middle2.heatFlow,
          tank1_VolumeCalculation2.heatFlow) annotation (Line(
            points={{-58.76,-8.9},{-58.76,-3.825},{-36.7273,-3.825}},
            color={0,0,0},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (Diagram(graphics),
          experiment(
            StopTime=10000,
            NumberOfIntervals=10000,
            Tolerance=1e-006),
          __Dymola_experimentSetupOutput);
      end TESTS_AREA_Inner;
    end TESTS;
    annotation (Icon(graphics={Bitmap(extent={{-80,-66},{90,64}}, fileName=
                "modelica://HySDeP/Graphics/PCM_tank.png")}), Diagram(graphics=
           {Bitmap(extent={{-74,-56},{80,84}}, fileName=
                "modelica://HySDeP/Graphics/PCM_tank.png")}));
  end CHG;
  annotation (uses(Modelica(version="3.2.1")),
                                             Icon(graphics={Bitmap(extent={{
              -102,-96},{92,96}}, fileName=
              "modelica://HySDeP/Graphics/Modeling_Platform.png")}),
    Diagram(graphics={Bitmap(extent={{-106,-86},{106,84}}, fileName=
              "modelica://HySDeP/Graphics/Modeling_Platform.png")}));
end HySDeP;
