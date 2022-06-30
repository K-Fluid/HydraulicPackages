package Accumulators
  model IdealAccumulator "Simple tank with inlet/outlet ports"
      import Modelica.Constants.pi;
  
    // Tank properties
    SI.Height level(stateSelect=StateSelect.prefer, start=level_start_eps)
        "Level height of tank";
    SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";
  
    // Tank geometry
    parameter SI.Height height "Height of tank";
    parameter SI.Area crossArea "Area of tank";
  
    // Ambient
    parameter Medium.AbsolutePressure p_ambient=system.p_ambient
        "Tank surface pressure"
      annotation(Dialog(tab = "Assumptions", group = "Ambient"));
    parameter Medium.Temperature T_ambient=system.T_ambient
        "Tank surface Temperature"
      annotation(Dialog(tab = "Assumptions", group = "Ambient"));
  
    // Initialization
    parameter SI.Height level_start(min=0) = 0.5*height
        "Start value of tank level"
      annotation(Dialog(tab="Initialization"));
  
    // Mass and energy balance, ports
    extends Modelica.Fluid.Vessels.BaseClasses.PartialLumpedVessel(
      final fluidVolume = V,
      final fluidLevel = level,
      final fluidLevel_max = height,
      final vesselArea = crossArea,
      heatTransfer(surfaceAreas={crossArea+2*sqrt(crossArea*pi)*level}),
      final initialize_p = false,
      final p_start = p_ambient);
  
    protected
    final parameter SI.Height level_start_eps = max(level_start, Modelica.Constants.eps);
  
  equation
    // Total quantities
    V = crossArea*level "Volume of fluid";
    medium.p = p_ambient;
  
    // Source termsEnergy balance
    if Medium.singleState or energyDynamics == Types.Dynamics.SteadyState then
      Wb_flow = 0
          "Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)";
    else
      Wb_flow = -p_ambient*der(V);
    end if;
  
    //Determine port properties
    for i in 1:nPorts loop
      vessel_ps_static[i] = max(0, level - portsData_height[i])*system.g*medium.d + p_ambient;
    end for;
  
  initial equation
    if massDynamics == Types.Dynamics.FixedInitial then
      level = level_start_eps;
    elseif massDynamics == Types.Dynamics.SteadyStateInitial then
      der(level) = 0;
    end if;
  
      annotation (defaultComponentName="tank",
        Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.2), graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={255,255,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.VerticalCylinder),
            Rectangle(
              extent=DynamicSelect({{-100,-100},{100,10}}, {{-100,-100},{100,(-100
                   + 200*level/height)}}),
              fillColor={85,170,255},
              fillPattern=FillPattern.VerticalCylinder),
            Line(points={{-100,100},{-100,-100},{100,-100},{100,100}}),
            Text(
              extent={{-95,60},{95,40}},
              textString="level ="),
            Text(
              extent={{-95,-24},{95,-44}},
              textString=DynamicSelect("%level_start", String(
                  level,
                  minimumLength=1,
                  significantDigits=2)))}),
        Documentation(info="<html>
  <p>
  Model of a tank that is open to the ambient at the fixed pressure
  <code>p_ambient</code>.
  </p>
  <p>
  The vector of connectors <strong>ports</strong> represents fluid ports at configurable heights, relative to the bottom of tank.
  Fluid can flow either out of or in to each port.
  </p>
  The following assumptions are made:
  <ul>
  <li>The tank is filled with a single or multiple-substance medium having a density higher than the density of the ambient medium.</li>
  <li>The fluid has uniform density, temperature and mass fractions</li>
  <li>No liquid is leaving the tank through the open top; the simulation breaks with an assertion if the liquid level growths over the height.</li>
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
  </html>", revisions="<html>
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
  </html>"));
  end IdealAccumulator;
end Accumulators;
