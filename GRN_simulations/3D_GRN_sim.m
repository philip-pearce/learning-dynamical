%create model
Mobj = sbiomodel('cell');
%add compartment
compObj = addcompartment(Mobj,'comp');
compObj.CapacityUnits = 'liter';
%Reactions
Robj1 = addreaction(Mobj,'A01 -> A01 + a');
Robj2 = addreaction(Mobj,'A00 -> A00 + a');
Robj3 = addreaction(Mobj,'A10 -> A10 + a');
Robj4 = addreaction(Mobj,'A11 -> A11 + a');
Robj5 = addreaction(Mobj,'B01 -> B01 + b');
Robj6 = addreaction(Mobj,'B00 -> B00 + b');
Robj7 = addreaction(Mobj,'B10 -> B10 + b');
Robj8 = addreaction(Mobj,'B11 -> B11 + b');
Robj9 = addreaction(Mobj,'C01 -> C01 + c');
Robj10 = addreaction(Mobj,'C00 -> C00 + c');
Robj11 = addreaction(Mobj,'C10 -> C10 + c');
Robj12 = addreaction(Mobj,'C11 -> C11 + c');

Robj13 = addreaction(Mobj,'a -> phi');
Robj14 = addreaction(Mobj,'b -> phi');
Robj15 = addreaction(Mobj,'c -> phi');

%Repression reactions
Robj16 = addreaction(Mobj,'A00 + 2 c <-> A01');
Robj17 = addreaction(Mobj,'A10 + 2 c <-> A11');
Robj18 = addreaction(Mobj,'B00 + 2 a <-> B01');
Robj19 = addreaction(Mobj,'B10 + 2 a <-> B11');
Robj20 = addreaction(Mobj,'C00 + 2 b <-> C01');
Robj21 = addreaction(Mobj,'C10 + 2 b <-> C11');

%Activation reactions
Robj22 = addreaction(Mobj,'A00 + 2 a <-> A10');
Robj23 = addreaction(Mobj,'A01 + 2 a <-> A11');
Robj24 = addreaction(Mobj,'B00 + 2 b <-> B10');
Robj25 = addreaction(Mobj,'B01 + 2 b <-> B11');
Robj26 = addreaction(Mobj,'C00 + 2 c <-> C10');
Robj27 = addreaction(Mobj,'C01 + 2 c <-> C11');

%initialise species amount
Mobj.Species(1).InitialAmount = 0;
Mobj.Species(2).InitialAmount = 0;
Mobj.Species(3).InitialAmount = 1;
Mobj.Species(4).InitialAmount = 0;
Mobj.Species(5).InitialAmount = 0;
Mobj.Species(6).InitialAmount = 0;
Mobj.Species(7).InitialAmount = 0;
Mobj.Species(8).InitialAmount = 1;
Mobj.Species(9).InitialAmount = 0;
Mobj.Species(10).InitialAmount = 0;
Mobj.Species(11).InitialAmount = 0;
Mobj.Species(12).InitialAmount = 0;
Mobj.Species(13).InitialAmount = 1;
Mobj.Species(14).InitialAmount = 0;
Mobj.Species(15).InitialAmount = 0;
Mobj.Species(16).InitialAmount = 0;

%Units of species amounts
for i = 1:16
Mobj.Species(i).InitialAmountUnits = 'molecule';
end
%Mass action kinetics
Kobj1 = addkineticlaw(Robj1,'MassAction');
Kobj2 = addkineticlaw(Robj2,'MassAction');
Kobj3 = addkineticlaw(Robj3,'MassAction');
Kobj4 = addkineticlaw(Robj4,'MassAction');
Kobj5 = addkineticlaw(Robj5,'MassAction');
Kobj6 = addkineticlaw(Robj6,'MassAction');
Kobj7 = addkineticlaw(Robj7,'MassAction');
Kobj8 = addkineticlaw(Robj8,'MassAction');
Kobj9 = addkineticlaw(Robj9,'MassAction');
Kobj10 = addkineticlaw(Robj10,'MassAction');
Kobj11 = addkineticlaw(Robj11,'MassAction');
Kobj12 = addkineticlaw(Robj12,'MassAction');
Kobj13 = addkineticlaw(Robj13,'MassAction');
Kobj14 = addkineticlaw(Robj14,'MassAction');
Kobj15 = addkineticlaw(Robj15,'MassAction');
Kobj16 = addkineticlaw(Robj16,'MassAction');
Kobj17 = addkineticlaw(Robj17,'MassAction');
Kobj18 = addkineticlaw(Robj18,'MassAction');
Kobj19 = addkineticlaw(Robj19,'MassAction');
Kobj20 = addkineticlaw(Robj20,'MassAction');
Kobj21 = addkineticlaw(Robj21,'MassAction');
Kobj22 = addkineticlaw(Robj22,'MassAction');
Kobj23 = addkineticlaw(Robj23,'MassAction');
Kobj24 = addkineticlaw(Robj24,'MassAction');
Kobj25 = addkineticlaw(Robj25,'MassAction');
Kobj26 = addkineticlaw(Robj26,'MassAction');
Kobj27 = addkineticlaw(Robj27,'MassAction');

%Rate parameters
g0 = 5;
g1 = 14;
k = 1;
hr = 1e-4;
fr = 1e-2;
ha = 2;
fa = 1e-1;
    %Rate parameter for Reaction 1
    Pobj1 = addparameter(Kobj1,'g0');
    Pobj1.Value = g0;
    Pobj1.ValueUnits = '1/second';
    Kobj1.ParameterVariableNames = 'g0';
    %Rate parameter for Reaction 2
    Pobj2 = addparameter(Kobj2,'g0');
    Pobj2.Value = g0;
    Pobj2.ValueUnits = '1/second';
    Kobj2.ParameterVariableNames = 'g0';
    %Rate parameter for Reaction 3
    Pobj3 = addparameter(Kobj3,'g1');
    Pobj3.Value = g1;
    Pobj3.ValueUnits = '1/second';
    Kobj3.ParameterVariableNames = 'g1';
    %Rate parameter for Reaction 4
    Pobj4 = addparameter(Kobj4,'g0');
    Pobj4.Value = g0;
    Pobj4.ValueUnits = '1/second';
    Kobj4.ParameterVariableNames = 'g0';
    %Rate parameter for Reaction 5
    Pobj5 = addparameter(Kobj5,'g0');
    Pobj5.Value = g0;
    Pobj5.ValueUnits = '1/second';
    Kobj5.ParameterVariableNames = 'g0';
    %Rate parameter for Reaction 6
    Pobj6 = addparameter(Kobj6,'g0');
    Pobj6.Value = g0;
    Pobj6.ValueUnits = '1/second';
    Kobj6.ParameterVariableNames = 'g0';
    %Rate parameter for Reaction 7
    Pobj7 = addparameter(Kobj7,'g1');
    Pobj7.Value = g1;
    Pobj7.ValueUnits = '1/second';
    Kobj7.ParameterVariableNames = 'g1';    
    %Rate parameter for Reaction 8
    Pobj8 = addparameter(Kobj8,'g0');
    Pobj8.Value = g0;
    Pobj8.ValueUnits = '1/second';
    Kobj8.ParameterVariableNames = 'g0';
    %Rate parameter for Reaction 9
    Pobj9 = addparameter(Kobj9,'g0');
    Pobj9.Value = g0;
    Pobj9.ValueUnits = '1/second';
    Kobj9.ParameterVariableNames = 'g0';
    %Rate parameter for Reaction 10
    Pobj10 = addparameter(Kobj10,'g0');
    Pobj10.Value = g0;
    Pobj10.ValueUnits = '1/second';
    Kobj10.ParameterVariableNames = 'g0';
    %Rate parameter for Reaction 11
    Pobj11 = addparameter(Kobj11,'g1');
    Pobj11.Value = g1;
    Pobj11.ValueUnits = '1/second';
    Kobj11.ParameterVariableNames = 'g1';    
    %Rate parameter for Reaction 12
    Pobj12 = addparameter(Kobj12,'g0');
    Pobj12.Value = g0;
    Pobj12.ValueUnits = '1/second';
    Kobj12.ParameterVariableNames = 'g0';
    
    %Rate parameter for Reaction 13
    Pobj13 = addparameter(Kobj13,'k');
    Pobj13.Value = k;
    Pobj13.ValueUnits = '1/second';
    Kobj13.ParameterVariableNames = 'k';    
    %Rate parameter for Reaction 14
    Pobj14 = addparameter(Kobj14,'k');
    Pobj14.Value = k;
    Pobj14.ValueUnits = '1/second';
    Kobj14.ParameterVariableNames = 'k';
    %Rate parameter for Reaction 15
    Pobj15 = addparameter(Kobj15,'k');
    Pobj15.Value = k;
    Pobj15.ValueUnits = '1/second';
    Kobj15.ParameterVariableNames = 'k';
    
    
    %Rate parameter for Reaction 16
    Pobj16 = addparameter(Kobj16,'hr');
    Pobj16.Value = hr;
    Pobj16.ValueUnits = '1/(molecule*molecule*second)';
    Pobj16r = addparameter(Kobj16,'fr');
    Pobj16r.Value = fr;
    Pobj16r.ValueUnits = '1/second';  
    Kobj16.ParameterVariableNames = {'hr','fr'};
    %Rate parameter for Reaction 17
    Pobj17 = addparameter(Kobj17,'hr');
    Pobj17.Value = hr;
    Pobj17.ValueUnits = '1/(molecule*molecule*second)';
    Pobj17r = addparameter(Kobj17,'fr');
    Pobj17r.Value = fr;
    Pobj17r.ValueUnits = '1/second';  
    Kobj17.ParameterVariableNames = {'hr','fr'};
    %Rate parameter for Reaction 18
    Pobj18 = addparameter(Kobj18,'hr');
    Pobj18.Value = hr;
    Pobj18.ValueUnits = '1/(molecule*molecule*second)';
    Pobj18r = addparameter(Kobj18,'fr');
    Pobj18r.Value = fr;
    Pobj18r.ValueUnits = '1/second';  
    Kobj18.ParameterVariableNames = {'hr','fr'};    
    %Rate parameter for Reaction 19
    Pobj19 = addparameter(Kobj19,'hr');
    Pobj19.Value = hr;
    Pobj19.ValueUnits = '1/(molecule*molecule*second)';
    Pobj19r = addparameter(Kobj19,'fr');
    Pobj19r.Value = fr;
    Pobj19r.ValueUnits = '1/second';  
    Kobj19.ParameterVariableNames = {'hr','fr'};    
    %Rate parameter for Reaction 20
    Pobj20 = addparameter(Kobj20,'hr');
    Pobj20.Value = hr;
    Pobj20.ValueUnits = '1/(molecule*molecule*second)';
    Pobj20r = addparameter(Kobj20,'fr');
    Pobj20r.Value = fr;
    Pobj20r.ValueUnits = '1/second';  
    Kobj20.ParameterVariableNames = {'hr','fr'};    
    %Rate parameter for Reaction 21
    Pobj21 = addparameter(Kobj21,'hr');
    Pobj21.Value = hr;
    Pobj21.ValueUnits = '1/(molecule*molecule*second)';
    Pobj21r = addparameter(Kobj21,'fr');
    Pobj21r.Value = fr;
    Pobj21r.ValueUnits = '1/second';  
    Kobj21.ParameterVariableNames = {'hr','fr'};   

    %Rate parameter for Reaction 22
    Pobj22 = addparameter(Kobj22,'ha');
    Pobj22.Value = ha;
    Pobj22.ValueUnits = '1/(molecule*molecule*second)';
    Pobj22r = addparameter(Kobj22,'fa');
    Pobj22r.Value = fa;
    Pobj22r.ValueUnits = '1/second';  
    Kobj22.ParameterVariableNames = {'ha','fa'};
    %Rate parameter for Reaction 23
    Pobj23 = addparameter(Kobj23,'ha');
    Pobj23.Value = ha;
    Pobj23.ValueUnits = '1/(molecule*molecule*second)';
    Pobj23r = addparameter(Kobj23,'fa');
    Pobj23r.Value = fa;
    Pobj23r.ValueUnits = '1/second';  
    Kobj23.ParameterVariableNames = {'ha','fa'};   
    %Rate parameter for Reaction 24
    Pobj24 = addparameter(Kobj24,'ha');
    Pobj24.Value = ha;
    Pobj24.ValueUnits = '1/(molecule*molecule*second)';
    Pobj24r = addparameter(Kobj24,'fa');
    Pobj24r.Value = fa;
    Pobj24r.ValueUnits = '1/second';  
    Kobj24.ParameterVariableNames = {'ha','fa'};   
    %Rate parameter for Reaction 25
    Pobj25 = addparameter(Kobj25,'ha');
    Pobj25.Value = ha;
    Pobj25.ValueUnits = '1/(molecule*molecule*second)';
    Pobj25r = addparameter(Kobj25,'fa');
    Pobj25r.Value = fa;
    Pobj25r.ValueUnits = '1/second';  
    Kobj25.ParameterVariableNames = {'ha','fa'};   
    %Rate parameter for Reaction 26
    Pobj26 = addparameter(Kobj26,'ha');
    Pobj26.Value = ha;
    Pobj26.ValueUnits = '1/(molecule*molecule*second)';
    Pobj26r = addparameter(Kobj26,'fa');
    Pobj26r.Value = fa;
    Pobj26r.ValueUnits = '1/second';  
    Kobj26.ParameterVariableNames = {'ha','fa'};   
    %Rate parameter for Reaction 27
    Pobj27 = addparameter(Kobj27,'ha');
    Pobj27.Value = ha;
    Pobj27.ValueUnits = '1/(molecule*molecule*second)';
    Pobj27r = addparameter(Kobj27,'fa');
    Pobj27r.Value = fa;
    Pobj27r.ValueUnits = '1/second';  
    Kobj27.ParameterVariableNames = {'ha','fa'};   
    
%simulate model stochastically
configset = getconfigset(Mobj);
configset.CompileOptions.DimensionalAnalysis = true;
configset.SolverType = 'ssa';
%configset.SolverType = 'ode15s';
configset.StopTime = 10000000;
solver = configset.SolverOptions;
%set(configset.SolverOptions, 'AbsoluteTolerance', 1.0e-8);
solver.LogDecimation = 100;

%do sim
[t, simdata, names] = sbiosimulate(Mobj);