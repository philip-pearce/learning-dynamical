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
Robj13 = addreaction(Mobj,'D01 -> D01 + d');
Robj14 = addreaction(Mobj,'D00 -> D00 + d');
Robj15 = addreaction(Mobj,'D10 -> D10 + d');
Robj16 = addreaction(Mobj,'D11 -> D11 + d');

%decay
Robj17 = addreaction(Mobj,'a -> phi');
Robj18 = addreaction(Mobj,'b -> phi');
Robj19 = addreaction(Mobj,'c -> phi');
Robj20 = addreaction(Mobj,'d -> phi');


%Repression reactions
Robj21 = addreaction(Mobj,'A00 + 2 d <-> A01');
Robj22 = addreaction(Mobj,'A10 + 2 d <-> A11');
Robj23 = addreaction(Mobj,'B00 + 2 a <-> B01');
Robj24 = addreaction(Mobj,'B10 + 2 a <-> B11');
Robj25 = addreaction(Mobj,'C00 + 2 b <-> C01');
Robj26 = addreaction(Mobj,'C10 + 2 b <-> C11');
Robj27 = addreaction(Mobj,'D00 + 2 c <-> D01');
Robj28 = addreaction(Mobj,'D10 + 2 c <-> D11');

%Activation reactions
Robj29 = addreaction(Mobj,'A00 + 2 a <-> A10');
Robj30 = addreaction(Mobj,'A01 + 2 a <-> A11');
Robj31 = addreaction(Mobj,'B00 + 2 b <-> B10');
Robj32 = addreaction(Mobj,'B01 + 2 b <-> B11');
Robj33 = addreaction(Mobj,'C00 + 2 c <-> C10');
Robj34 = addreaction(Mobj,'C01 + 2 c <-> C11');
Robj35 = addreaction(Mobj,'D00 + 2 d <-> D10');
Robj36 = addreaction(Mobj,'D01 + 2 d <-> D11');

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
Mobj.Species(17).InitialAmount = 0;
Mobj.Species(18).InitialAmount = 1;
Mobj.Species(19).InitialAmount = 0;
Mobj.Species(20).InitialAmount = 0;
Mobj.Species(21).InitialAmount = 0;

%Units of species amounts
for i = 1:21
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
Kobj28 = addkineticlaw(Robj28,'MassAction');
Kobj29 = addkineticlaw(Robj29,'MassAction');
Kobj30 = addkineticlaw(Robj30,'MassAction');
Kobj31 = addkineticlaw(Robj31,'MassAction');
Kobj32 = addkineticlaw(Robj32,'MassAction');
Kobj33 = addkineticlaw(Robj33,'MassAction');
Kobj34 = addkineticlaw(Robj34,'MassAction');
Kobj35 = addkineticlaw(Robj35,'MassAction');
Kobj36 = addkineticlaw(Robj36,'MassAction');

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
    Pobj13 = addparameter(Kobj13,'g0');
    Pobj13.Value = g0;
    Pobj13.ValueUnits = '1/second';
    Kobj13.ParameterVariableNames = 'g0';
    %Rate parameter for Reaction 14
    Pobj14 = addparameter(Kobj14,'g0');
    Pobj14.Value = g0;
    Pobj14.ValueUnits = '1/second';
    Kobj14.ParameterVariableNames = 'g0';
    %Rate parameter for Reaction 15
    Pobj15 = addparameter(Kobj15,'g1');
    Pobj15.Value = g1;
    Pobj15.ValueUnits = '1/second';
    Kobj15.ParameterVariableNames = 'g1';    
    %Rate parameter for Reaction 16
    Pobj16 = addparameter(Kobj16,'g0');
    Pobj16.Value = g0;
    Pobj16.ValueUnits = '1/second';
    Kobj16.ParameterVariableNames = 'g0';
    
    
    
    %Rate parameter for Reaction 17
    Pobj17 = addparameter(Kobj17,'k');
    Pobj17.Value = k;
    Pobj17.ValueUnits = '1/second';
    Kobj17.ParameterVariableNames = 'k';    
    %Rate parameter for Reaction 18
    Pobj18 = addparameter(Kobj18,'k');
    Pobj18.Value = k;
    Pobj18.ValueUnits = '1/second';
    Kobj18.ParameterVariableNames = 'k';
    %Rate parameter for Reaction 19
    Pobj19 = addparameter(Kobj19,'k');
    Pobj19.Value = k;
    Pobj19.ValueUnits = '1/second';
    Kobj19.ParameterVariableNames = 'k';
    %Rate parameter for Reaction 20
    Pobj20 = addparameter(Kobj20,'k');
    Pobj20.Value = k;
    Pobj20.ValueUnits = '1/second';
    Kobj20.ParameterVariableNames = 'k';
    
    
    %Rate parameter for Reaction 21
    Pobj21 = addparameter(Kobj21,'hr');
    Pobj21.Value = hr;
    Pobj21.ValueUnits = '1/(molecule*molecule*second)';
    Pobj21r = addparameter(Kobj21,'fr');
    Pobj21r.Value = fr;
    Pobj21r.ValueUnits = '1/second';  
    Kobj21.ParameterVariableNames = {'hr','fr'};
    %Rate parameter for Reaction 22
    Pobj22 = addparameter(Kobj22,'hr');
    Pobj22.Value = hr;
    Pobj22.ValueUnits = '1/(molecule*molecule*second)';
    Pobj22r = addparameter(Kobj22,'fr');
    Pobj22r.Value = fr;
    Pobj22r.ValueUnits = '1/second';  
    Kobj22.ParameterVariableNames = {'hr','fr'};
    %Rate parameter for Reaction 23
    Pobj23 = addparameter(Kobj23,'hr');
    Pobj23.Value = hr;
    Pobj23.ValueUnits = '1/(molecule*molecule*second)';
    Pobj23r = addparameter(Kobj23,'fr');
    Pobj23r.Value = fr;
    Pobj23r.ValueUnits = '1/second';  
    Kobj23.ParameterVariableNames = {'hr','fr'};
    %Rate parameter for Reaction 24
    Pobj24 = addparameter(Kobj24,'hr');
    Pobj24.Value = hr;
    Pobj24.ValueUnits = '1/(molecule*molecule*second)';
    Pobj24r = addparameter(Kobj24,'fr');
    Pobj24r.Value = fr;
    Pobj24r.ValueUnits = '1/second';  
    Kobj24.ParameterVariableNames = {'hr','fr'};
    %Rate parameter for Reaction 25
    Pobj25 = addparameter(Kobj25,'hr');
    Pobj25.Value = hr;
    Pobj25.ValueUnits = '1/(molecule*molecule*second)';
    Pobj25r = addparameter(Kobj25,'fr');
    Pobj25r.Value = fr;
    Pobj25r.ValueUnits = '1/second';  
    Kobj25.ParameterVariableNames = {'hr','fr'};
    %Rate parameter for Reaction 26
    Pobj26 = addparameter(Kobj26,'hr');
    Pobj26.Value = hr;
    Pobj26.ValueUnits = '1/(molecule*molecule*second)';
    Pobj26r = addparameter(Kobj26,'fr');
    Pobj26r.Value = fr;
    Pobj26r.ValueUnits = '1/second';  
    Kobj26.ParameterVariableNames = {'hr','fr'};
    %Rate parameter for Reaction 27
    Pobj27 = addparameter(Kobj27,'hr');
    Pobj27.Value = hr;
    Pobj27.ValueUnits = '1/(molecule*molecule*second)';
    Pobj27r = addparameter(Kobj27,'fr');
    Pobj27r.Value = fr;
    Pobj27r.ValueUnits = '1/second';  
    Kobj27.ParameterVariableNames = {'hr','fr'};
    %Rate parameter for Reaction 28
    Pobj28 = addparameter(Kobj28,'hr');
    Pobj28.Value = hr;
    Pobj28.ValueUnits = '1/(molecule*molecule*second)';
    Pobj28r = addparameter(Kobj28,'fr');
    Pobj28r.Value = fr;
    Pobj28r.ValueUnits = '1/second';  
    Kobj28.ParameterVariableNames = {'hr','fr'};

    %Rate parameter for Reaction 29
    Pobj29 = addparameter(Kobj29,'ha');
    Pobj29.Value = ha;
    Pobj29.ValueUnits = '1/(molecule*molecule*second)';
    Pobj29r = addparameter(Kobj29,'fa');
    Pobj29r.Value = fa;
    Pobj29r.ValueUnits = '1/second';  
    Kobj29.ParameterVariableNames = {'ha','fa'};
    %Rate parameter for Reaction 30
    Pobj30 = addparameter(Kobj30,'ha');
    Pobj30.Value = ha;
    Pobj30.ValueUnits = '1/(molecule*molecule*second)';
    Pobj30r = addparameter(Kobj30,'fa');
    Pobj30r.Value = fa;
    Pobj30r.ValueUnits = '1/second';  
    Kobj30.ParameterVariableNames = {'ha','fa'};   
    %Rate parameter for Reaction 31
    Pobj31 = addparameter(Kobj31,'ha');
    Pobj31.Value = ha;
    Pobj31.ValueUnits = '1/(molecule*molecule*second)';
    Pobj31r = addparameter(Kobj31,'fa');
    Pobj31r.Value = fa;
    Pobj31r.ValueUnits = '1/second';  
    Kobj31.ParameterVariableNames = {'ha','fa'};   
    %Rate parameter for Reaction 32
    Pobj32 = addparameter(Kobj32,'ha');
    Pobj32.Value = ha;
    Pobj32.ValueUnits = '1/(molecule*molecule*second)';
    Pobj32r = addparameter(Kobj32,'fa');
    Pobj32r.Value = fa;
    Pobj32r.ValueUnits = '1/second';  
    Kobj32.ParameterVariableNames = {'ha','fa'};   
    %Rate parameter for Reaction 33
    Pobj33 = addparameter(Kobj33,'ha');
    Pobj33.Value = ha;
    Pobj33.ValueUnits = '1/(molecule*molecule*second)';
    Pobj33r = addparameter(Kobj33,'fa');
    Pobj33r.Value = fa;
    Pobj33r.ValueUnits = '1/second';  
    Kobj33.ParameterVariableNames = {'ha','fa'};   
    %Rate parameter for Reaction 34
    Pobj34 = addparameter(Kobj34,'ha');
    Pobj34.Value = ha;
    Pobj34.ValueUnits = '1/(molecule*molecule*second)';
    Pobj34r = addparameter(Kobj34,'fa');
    Pobj34r.Value = fa;
    Pobj34r.ValueUnits = '1/second';  
    Kobj34.ParameterVariableNames = {'ha','fa'};   
    %Rate parameter for Reaction 35
    Pobj35 = addparameter(Kobj35,'ha');
    Pobj35.Value = ha;
    Pobj35.ValueUnits = '1/(molecule*molecule*second)';
    Pobj35r = addparameter(Kobj35,'fa');
    Pobj35r.Value = fa;
    Pobj35r.ValueUnits = '1/second';  
    Kobj35.ParameterVariableNames = {'ha','fa'};   
    %Rate parameter for Reaction 36
    Pobj36 = addparameter(Kobj36,'ha');
    Pobj36.Value = ha;
    Pobj36.ValueUnits = '1/(molecule*molecule*second)';
    Pobj36r = addparameter(Kobj36,'fa');
    Pobj36r.Value = fa;
    Pobj36r.ValueUnits = '1/second';  
    Kobj36.ParameterVariableNames = {'ha','fa'};   
    
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