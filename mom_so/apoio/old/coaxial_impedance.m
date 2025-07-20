function [R, L] = coaxial_impedance(frequencies, a, b, c, sigma)
    % Constants
    MU0 = 4 * pi * 1e-7; % Permeability of free space
    
    % Preallocate the arrays for inductance and resistance
    L = zeros(size(frequencies));
    R = zeros(size(frequencies));

    for f = 1:length(frequencies)
        w = 2 * pi * frequencies(f);
        jw = 1j * w;
        
        % Calculate L'
        L_prime = MU0 / (2 * pi) * log(b / a);
        
        % Calculate eta and gamma_c
        eta = sqrt(jw * MU0 / sigma);
        gamma_c = sqrt(jw * MU0 * sigma);
        
        % Calculate Za(omega)
        I0_gamma_a = besseli(0, gamma_c * a);
        I1_gamma_a = besseli(1, gamma_c * a);
        Za = eta / (2 * pi * a) * (I0_gamma_a / I1_gamma_a);
        
        % Calculate Zb(omega) using the revised formula
        I0_gamma_c_b = besseli(0, gamma_c * b);
        I1_gamma_c_b = besseli(1, gamma_c * b);
        K0_gamma_c_b = besselk(0, gamma_c * b);
        K1_gamma_c_b = besselk(1, gamma_c * b);
        I1_gamma_c_c = besseli(1, gamma_c * c);
        K1_gamma_c_c = besselk(1, gamma_c * c);
        
        Zb = (eta / (2 * pi * b)) * ...
            ((I0_gamma_c_b * K1_gamma_c_c + K0_gamma_c_b * I1_gamma_c_c) / ...
            (I1_gamma_c_c * K1_gamma_c_b - I1_gamma_c_b * K1_gamma_c_c));
        
        % Total impedance
        Z = jw * L_prime + Za + Zb;

        % Separate real and imaginary parts
        L(f) = imag(Z) / w;       % Inductance
        R(f) = real(Z);           % Resistance
    end
end


