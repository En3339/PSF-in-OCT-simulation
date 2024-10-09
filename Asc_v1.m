function result = Asc(Mat_Sc, zr, kvec, Ai)
    % Constants
    W0 = 100e-6;  % Waist radius
    Nx = 32;
    dx = 0.2e-6;
    dy = dx;
    
    % Other constants and parameters
    dx_sc = 1e-7;
    dz_sc = 5e-7;
    numelzs = sum(sum(Mat_Sc == 200));
    
    % ... (other variable initializations)
    
    Esc = zeros(Nx, Nx, numel(kvec)); % scattered field
    Er = zeros(Nx, Nx, numel(kvec));  % reference field
    
    % ... (rest of your code)
    
    % Plot the result
    A = fft(Ik);
    zbar = (0:(numel(A)-1)) / numel(A) / diff(kvec(10:11)/2/pi);
    plot(zbar, abs(A));
    
    % Return the result if needed
    result = A;
end
