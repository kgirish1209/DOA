clear all;
close all;

display ("DOA Estimation Analysis for Conventional/Barlett and Capon Beam forming techniques");
display("Assumes two Sources incident at an angle [theta1_degree, theta2_degree");
display("Provide inputs M = # of antenna elements on ULA and [], true incident angles");

doa_Algo_fn(2,[-15 25]);

doa_Algo_fn(4,[-15 25]);

doa_Algo_fn(8,[-15 25]);
doa_Algo_fn(16,[-15 25]);
doa_Algo_fn(8,[-5 5]);
hold on;
