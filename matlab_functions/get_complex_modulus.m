function [Gprime, Gprimeprime,omega1,omega2] =get_complex_modulus(stress,strain,omega,dt)
    stress = stress(round(length(stress)/2):end);
    strain = strain(round(length(strain)/2):end);  
    
    sampling_freq = 1/dt;
    frequency = sampling_freq*(0:length(strain)-1)'/length(strain);
    sigma = fft(stress);
    sigma = sigma';
    gamma = fft(strain);
    gamma=gamma';

    [~,k]=maxk(abs(gamma(1:end)),3);
    [~,id] = min(abs(frequency(k)*2*pi-omega));
    k1=k(id);
    omega1=frequency(k1)*2*pi;

    [~,k1]=min(abs(frequency*2*pi-omega));
    omega2=frequency(k1)*2*pi;

    phi_g =angle(gamma(k1));
    phi_s = angle(sigma(k1));
    delta = phi_g-phi_s;

    Gprime = abs(sigma(k1))/abs(gamma(k1))*cos(delta);
    Gprimeprime = abs(sigma(k1))/abs(gamma(k1))*sin(delta);
    
  
end