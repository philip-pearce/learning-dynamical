function [mfpts,statdist] = MFPTmatrixWithHessianSCALED(state_energies, state_hessians, saddle_energies, saddle_hessians, transitions, temperature, friction)
%MFPTMATRIX Returns :
%           mfpts - the matrix such that (i,j) is the MFPT from i to j
%           statdist - the stationary distribution of the Markov chain from
%           the Kramers rates
% Inputs:   
%           state_energies - the energy of each state
%           state_hessians - entry (m,:,:) is the Hessian of the energy at
%           minimum m
%           saddle_energies - the energy of each saddle
%           saddle_hessians - entry (s,:,:) is the Hessian of the energy at
%           saddle s
%           transitions - a 3-column matrix where each row is
%                            [state1 saddle state2]
%                         for each state transition pathway found
%           temperature - the temperature of the thermal bath
%           friction - the damping of the diffusive walk

N_states = length(state_energies);
N_saddles = length(saddle_energies);
N_dims = size(saddle_hessians,2);

% compute eigenvalues of the Hessians
state_eigenvalues = zeros(N_states, N_dims);
saddle_eigenvalues = zeros(N_saddles, N_dims);
for i = 1:N_states
    state_eigenvalues(i,:) = eig(squeeze(state_hessians(i,:,:)));
end
for i = 1:N_saddles
    saddle_eigenvalues(i,:) = eig(squeeze(saddle_hessians(i,:,:)));
end

% precompute the elements necessary to form the prefactors
% stable states: product of angular frequencies (sqrt of |determinant| of Hessian)
% [n.b.  omega^2 =  -eigenvalue ]
state_prefactor = sqrt(abs(prod(state_eigenvalues,2)));

% saddles: unstable angular frequency / product of stable angular frequencies
saddle_unstable_evals = zeros(N_saddles,1);
saddle_stable_evals = zeros(N_saddles,N_dims-1);
for i = 1:N_saddles
    [~,idxs] = min(saddle_eigenvalues(i,:));
%     if length(idxs) ~= 1
%         error('Saddle with <>1 unstable eigenvalue');
%     end
    saddle_unstable_evals(i) = saddle_eigenvalues(i,idxs);
    %quick way to avoid errors - if you have any small negative
    %eigenvalues just make them positive
    saddle_stable_evals(i,:) = abs(saddle_eigenvalues(i,(1:N_dims) ~= idxs));
end
saddle_prefactor = sqrt(abs(saddle_unstable_evals ./ prod(saddle_stable_evals,2)));

% form the Markov matrix P
P = zeros(N_states, N_states);
for i = 1:size(transitions,1)
    state1 = transitions(i,1);
    saddle = transitions(i,2);
    state2 = transitions(i,3);
    
    % transitions happen at rate prefactor*exp(-EnergyBarrier / Temperature)
    prefactor = state_prefactor(state1)*saddle_prefactor(saddle) / (2*pi*friction);
    deltaE12 = saddle_energies(saddle) - state_energies(state1);
    assert(deltaE12 > 0);
    P(state1,state2) = prefactor*exp(-deltaE12 / temperature);
    
    % ... and in the opposite direction too
    prefactor = state_prefactor(state2)*saddle_prefactor(saddle) / (2*pi*friction);
    deltaE21 = saddle_energies(saddle) - state_energies(state2);
    assert(deltaE21 > 0);
    P(state2,state1) = prefactor*exp(-deltaE21 / temperature);
end

% the rate matrix must have its rows sum to zero
P = P - diag(sum(P,2));

% get the MFPTs from every vertex to each vertex as a target
mfpts = zeros(N_states,N_states);
for i = 1:N_states
    mfpts(:,i) = TargetMFPT(P, i);
end

% compute the stationary distribution v arising from P
% v satisfies v.P = 0 (i.e. is the left eigenvector for eigenvalue 0)
% all other eigenvalues are < 0
[V,~] = eigs(P',1,'lr');
statdist = V / sum(V);


function [ mfpt_target ] = TargetMFPT( P, target_state )
%TARGETMFPT Compute the mean first passage times from each state to the
%target state with transition matrix P

% set up the linear system whose solution vector k_i is the MFPT to target
% when started at i
% (see Norris, Thm 1.3.5)
LHS = -P;
LHS(target_state,:) = 0;
LHS(target_state,target_state) = 1;
RHS = ones(N_states,1);
RHS(target_state) = 0;

% solve it for the MFPTs
mfpt_target = LHS \ RHS;

end

end