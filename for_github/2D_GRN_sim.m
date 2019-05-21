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
Robj9 = addreaction(Mobj,'a -> phi');
Robj10 = addreaction(Mobj,'b -> phi');
Robj11 = addreaction(Mobj,'A00 + 2 b <-> A01');
Robj12 = addreaction(Mobj,'A10 + 2 b <-> A11');
Robj13 = addreaction(Mobj,'B00 + 2 a <-> B01');
Robj14 = addreaction(Mobj,'B10 + 2 a <-> B11');
Robj15 = addreaction(Mobj,'A00 + 2 a <-> A10');
Robj16 = addreaction(Mobj,'A01 + 2 a <-> A11');
Robj17 = addreaction(Mobj,'B00 + 2 b <-> B10');
Robj18 = addreaction(Mobj,'B01 + 2 b <-> B11');
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
%Units of species amounts
for i = 1:11
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
%Rate parameters
g0 = 4;
g1 = 16;
k=1;
hr = 1e-4;
fr = 1e-2;
ha = 1e-1;
fa = 1;
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
    Pobj9 = addparameter(Kobj9,'k');
    Pobj9.Value = k;
    Pobj9.ValueUnits = '1/second';
    Kobj9.ParameterVariableNames = 'k';
    %Rate parameter for Reaction 10
    Pobj10 = addparameter(Kobj10,'k');
    Pobj10.Value = k;
    Pobj10.ValueUnits = '1/second';
    Kobj10.ParameterVariableNames = 'k';    
    %Rate parameter for Reaction 11
    Pobj11 = addparameter(Kobj11,'hr');
    Pobj11.Value = hr;
    Pobj11.ValueUnits = '1/(molecule*molecule*second)';
    Pobj11r = addparameter(Kobj11,'fr');
    Pobj11r.Value = fr;
    Pobj11r.ValueUnits = '1/second';  
    Kobj11.ParameterVariableNames = {'hr','fr'};
    %Rate parameter for Reaction 12
    Pobj12 = addparameter(Kobj12,'hr');
    Pobj12.Value = hr;
    Pobj12.ValueUnits = '1/(molecule*molecule*second)';
    Pobj12r = addparameter(Kobj12,'fr');
    Pobj12r.Value = fr;
    Pobj12r.ValueUnits = '1/second';  
    Kobj12.ParameterVariableNames = {'hr','fr'};   
    %Rate parameter for Reaction 13
    Pobj13 = addparameter(Kobj13,'hr');
    Pobj13.Value = hr;
    Pobj13.ValueUnits = '1/(molecule*molecule*second)';
    Pobj13r = addparameter(Kobj13,'fr');
    Pobj13r.Value = fr;
    Pobj13r.ValueUnits = '1/second';  
    Kobj13.ParameterVariableNames = {'hr','fr'};
    %Rate parameter for Reaction 14
    Pobj14 = addparameter(Kobj14,'hr');
    Pobj14.Value = hr;
    Pobj14.ValueUnits = '1/(molecule*molecule*second)';
    Pobj14r = addparameter(Kobj14,'fr');
    Pobj14r.Value = fr;
    Pobj14r.ValueUnits = '1/second';  
    Kobj14.ParameterVariableNames = {'hr','fr'};   
    %Rate parameter for Reaction 15
    Pobj15 = addparameter(Kobj15,'ha');
    Pobj15.Value = ha;
    Pobj15.ValueUnits = '1/(molecule*molecule*second)';
    Pobj15r = addparameter(Kobj15,'fa');
    Pobj15r.Value = fa;
    Pobj15r.ValueUnits = '1/second';  
    Kobj15.ParameterVariableNames = {'ha','fa'};
    %Rate parameter for Reaction 16
    Pobj16 = addparameter(Kobj16,'ha');
    Pobj16.Value = ha;
    Pobj16.ValueUnits = '1/(molecule*molecule*second)';
    Pobj16r = addparameter(Kobj16,'fa');
    Pobj16r.Value = fa;
    Pobj16r.ValueUnits = '1/second';  
    Kobj16.ParameterVariableNames = {'ha','fa'};   
    %Rate parameter for Reaction 17
    Pobj17 = addparameter(Kobj17,'ha');
    Pobj17.Value = ha;
    Pobj17.ValueUnits = '1/(molecule*molecule*second)';
    Pobj17r = addparameter(Kobj17,'fa');
    Pobj17r.Value = fa;
    Pobj17r.ValueUnits = '1/second';  
    Kobj17.ParameterVariableNames = {'ha','fa'};
    %Rate parameter for Reaction 18
    Pobj18 = addparameter(Kobj18,'ha');
    Pobj18.Value = ha;
    Pobj18.ValueUnits = '1/(molecule*molecule*second)';
    Pobj18r = addparameter(Kobj18,'fa');
    Pobj18r.Value = fa;
    Pobj18r.ValueUnits = '1/second';  
    Kobj18.ParameterVariableNames = {'ha','fa'};    

    
%simulate model stochastically
configset = getconfigset(Mobj);
configset.CompileOptions.DimensionalAnalysis = true;
configset.SolverType = 'ssa';
configset.StopTime = 1000000;
solver = configset.SolverOptions;
solver.LogDecimation = 100;
%configset.CompileOptions.DimensionalAnalysis = false;


%do sim
[t, simdata, names] = sbiosimulate(Mobj);