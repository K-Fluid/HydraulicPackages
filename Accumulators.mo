package Accumulators
  model IdealAccumulator "Accumulator model filled with ideal Nitrogen gas"
    import Modelica.Constants.pi;
    // スコープの調整
    import SI = Modelica.Units.SI;
    import Modelica.Fluid.Types;
    import Modelica.Fluid;
    // 封入気体の物性モデルと状態変数
    // Filled gas
    replaceable package N2 = Modelica.Media.IdealGases.SingleGases.N2;
    N2.BaseProperties gas;
    SI.Volume V_gas;
    SI.Mass M_gas;
    SI.Energy U_gas;
    // Accumulator properties
    SI.Volume V_oil(stateSelect = StateSelect.prefer, start = V_oil_start_eps) "Oil volume";
    // Accumulator geometry
    parameter SI.Volume V_total(min = pi / 6 * id ^ 3) = 0.001 "Volume of Shell";
    parameter SI.Diameter id = 0.1 "Inner Diameter of Shell";
    parameter SI.Length L_diaphragm = 0.1 "Length of Shell separated by Diaphragm";
    // Ambient
    parameter Medium.AbsolutePressure p_ambient = system.p_ambient "Ambient Pressure" annotation(
      Dialog(tab = "Assumptions", group = "Ambient"));
    parameter Medium.Temperature T_ambient = system.T_ambient "Ambient Temperature" annotation(
      Dialog(tab = "Assumptions", group = "Ambient"));
    // Initialization
    // 初期の油量
    parameter SI.Volume V_oil_start(min = Modelica.Constants.eps) = 1e-06 "Start value of oil volume" annotation(
      Dialog(tab = "Initialization"));
    parameter SI.Pressure p_filled = 5e+06 "Filled gas pressure";
    // Mass and energy balance, ports
    extends Modelica.Fluid.Vessels.BaseClasses.PartialLumpedVessel(
      final fluidVolume = V_oil, 
      heatTransfer(surfaceAreas = {pi / 2 * id ^ 2 + pi / 4 * id ^ 2 * (L_diaphragm - id / 2)}), 
      final initialize_p = false, 
      energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, 
      massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial,
      nPorts = 2, use_portsData = false);
    //    portsData = {
    //      Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.010, height = 0),
    //      Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.010, height = 0)},
  protected
    final parameter SI.Volume V_oil_start_eps = max(V_oil_start, Modelica.Constants.eps);
  equation
// 境界条件を追加
    if medium.p > p_filled then
      medium.p = gas.p;
    else
      gas.p = p_filled;
    end if;
// 空気の質量と内部エネルギーを計算
    V_gas = V_total - V_oil;
    M_gas = gas.d * V_gas;
    U_gas = gas.u * M_gas;
//　空気の質量保存則とエネルギー保存則
    der(M_gas) = 0;
// Source termsEnergy balance
    if Medium.singleState or energyDynamics == Types.Dynamics.SteadyState then
      Wb_flow = 0 "Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if oil volume changes)";
      der(U_gas) = 0;
    else
      Wb_flow = -medium.p * der(V_oil);
      der(U_gas) = -Wb_flow;
    end if;
//Determine port properties
// 静圧の計算
    for i in 1:nPorts loop
      vessel_ps_static[i] = medium.p;
    end for;
  initial equation
    if massDynamics == Types.Dynamics.FixedInitial then
      V_oil = V_oil_start_eps;
    elseif massDynamics == Types.Dynamics.SteadyStateInitial then
      der(V_oil) = 0;
    end if;
    gas.p = p_filled;
    medium.p = p_ambient;
    gas.T = T_ambient;
    annotation(
      Dialog(tab = "Initialization"));
    annotation(
      defaultComponentName = "Accumulator",
      Icon(graphics = {Rectangle(lineThickness = 2, extent = {{-50, 100}, {50, -100}}, radius = 50), Line(points = {{-50, 0}, {50, 0}}, thickness = 2), Polygon(lineThickness = 2, points = {{11.42, 19.78}, {0, 0}, {-11.42, 19.78}, {11.42, 19.78}})}),
      Documentation(info = "<html><head></head><body><p>
  Model of an accumulator that is filled with ideal nitrogen gas at the pressure
  <code>p_filled</code>.
  </p>
  <p>
  The vector of connectors <strong>ports</strong> represents fluid ports of the accumulator.
  Fluid can flow either out of or in to each port.
  </p>
  The following assumptions are made:
  <ul>
  <li>The accumulator is filled with a single or multiple-substance medium.</li>
  <li>The fluid has uniform density, temperature and mass fractions.</li>
  </ul>
  <p>
  The port pressures represent the pressures just after the outlet (or just before the inlet) in the attached pipe.
  The hydraulic resistances <code>portsData.zeta_in</code> and <code>portsData.zeta_out</code> determine the dissipative pressure drop between tank and port depending on
  the direction of mass flow. See <a href=\"modelica://Modelica.Fluid.Vessels.BaseClasses.VesselPortsData\">VesselPortsData</a> and <em>[Idelchik, Handbook of Hydraulic Resistance, 2004]</em>.
  </p>
  <p>
  With the setting <code>use_portsData=false</code>, the port pressure represents the static head
  at the height of the respective port.
  The relationship between pressure drop and mass flow rate at the port must then be provided by connected components;
  Heights of ports as well as kinetic and potential energy of fluid entering or leaving are not taken into account anymore.
  </p>
  <img src=\"modelica://Accumulators/Resources/Images/IdealAccumulator_1.png\">
  <img src=\"modelica://Accumulators/Resources/Images/IdealAccumulator_2.png\">
  </body></html>", revisions = "<html>
  <ul>
  <li><em>Dec. 12, 2008</em> by R&uuml;diger Franke: move port definitions
     to BaseClasses.PartialLumpedVessel; also use energy and mass balance from common base class</li>
  <li><em>Dec. 8, 2008</em> by Michael Wetter (LBNL):<br>
  Implemented trace substances.</li>
  <li><em>Jan. 6, 2006</em> by Katja Poschlad, Manuel Remelhe (AST Uni Dortmund),
     Martin Otter (DLR):<br>
     Implementation based on former tank model.</li>
  <li><em>Oct. 29, 2007</em> by Carsten Heinrich (ILK Dresden):<br>
  Adapted to the new fluid library interfaces:
  <ul> <li>FluidPorts_b is used instead of FluidPort_b (due to it is defined as an array of ports)</li>
      <li>Port name changed from port to ports</li></ul>Updated documentation.</li>
  <li><em>Apr. 25, 2006</em> by Katrin Pr&ouml;l&szlig; (TUHH):<br>
  Limitation to bottom ports only, added inlet and outlet loss factors.</li>
  </ul>
  </html>", __OpenModelica_infoHeader = "<html><head></head><body><h4><br></h4></body></html>"));
  end IdealAccumulator;

  model AccumulatorTest1
    replaceable package Medium = Modelica.Media.Water.StandardWater;
    Modelica.Blocks.Sources.Ramp ramp1(duration = 10, height = 0.1, offset = 0, startTime = 0) annotation(
      Placement(visible = true, transformation(origin = {-70, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {50, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp2(duration = 10, height = -0.2, offset = 0, startTime = 6) annotation(
      Placement(visible = true, transformation(origin = {60, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.MassFlowSource_T boundary(redeclare package Medium = Medium, nPorts = 1, use_m_flow_in = true) annotation(
      Placement(visible = true, transformation(origin = {-32, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.MassFlowSource_T boundary1(redeclare package Medium = Medium, m_flow = 0, nPorts = 1, use_m_flow_in = true) annotation(
      Placement(visible = true, transformation(origin = {30, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    IdealAccumulator acc(redeclare package Medium = Medium, V_oil(fixed = false), nPorts = 2) annotation(
      Placement(visible = true, transformation(origin = {0, 2}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  equation
    connect(ramp1.y, boundary.m_flow_in) annotation(
      Line(points = {{-58, 0}, {-52, 0}, {-52, -42}, {-42, -42}}, color = {0, 0, 127}));
    connect(ramp2.y, boundary1.m_flow_in) annotation(
      Line(points = {{72, 0}, {82, 0}, {82, -46}, {42, -46}}, color = {0, 0, 127}));
    connect(boundary.ports[1], acc.ports[1]) annotation(
      Line(points = {{-22, -50}, {0, -50}, {0, -18}}, color = {0, 127, 255}));
    connect(boundary1.ports[1], acc.ports[2]) annotation(
      Line(points = {{20, -50}, {0, -50}, {0, -18}}, color = {0, 127, 255}));
  end AccumulatorTest1;

  model AccumulatorTest2
    replaceable package Medium = Modelica.Media.Water.StandardWater;
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {50, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.MassFlowSource_T boundary(redeclare package Medium = Medium, nPorts = 1, use_m_flow_in = true) annotation(
      Placement(visible = true, transformation(origin = {-32, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    IdealAccumulator acc(redeclare package Medium = Medium, V_oil(fixed = false), nPorts = 1) annotation(
      Placement(visible = true, transformation(origin = {0, 2}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Sources.TimeTable timeTable(table = [0, 0; 4, .05; 5, .05; 12, -.05; 13, -.05; 17, 0])  annotation(
      Placement(visible = true, transformation(origin = {-80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(boundary.ports[1], acc.ports[1]) annotation(
      Line(points = {{-22, -50}, {0, -50}, {0, -18}}, color = {0, 127, 255}));
    connect(timeTable.y, boundary.m_flow_in) annotation(
      Line(points = {{-68, 0}, {-60, 0}, {-60, -46}, {-44, -46}}, color = {0, 0, 127}));
  end AccumulatorTest2;
  annotation(
    uses(Modelica(version = "4.0.0")));
end Accumulators;
