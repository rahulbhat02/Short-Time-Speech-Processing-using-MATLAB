% This function is used for pre-processing of the sound files
% pre-processing involves following steps
% 
% 1. Windowing
% 2. Pre-emphasis

function [p_data] = pre_process(x)

% x : input signal (could be in rows of frames i.e. row1 = frame1)

% applying window on each frame
[r,c] = size(x);

for i = 1 : r
x1(i,:) = x(i,:).*hamming(length(x(i,:)))'; % hamming window

end

% Apply a pre-emphasis filter. The pre-emphasis filter is a highpass
preemph = [1 -0.95];
x1 = filter(preemph, 1, x1);

% final output
p_data = x1;
end