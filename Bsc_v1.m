%zs = 200e-6 + rand(1,2)*10^(-4);
zs = 200e-6 + [100e-6 ];
zr = 0;
E0=1;
W_0 = 15e-6;           % Waist radius
kvec = linspace( 2*pi/(1300e-9+150e-9), 2*pi/(1300e-9-150e-9), 2048);%vector of wave numbers
E_sc = zeros(size(kvec));%scattered field
Er = zeros(size(kvec));%reference field

var_dx = 1e-6;
var_x = (-10:10)*var_dx;




B = zeros(20,2048);
for j = 1:length(var_x)
    Nx = 256;               
    dx = 0.5e-6;dy=dx;
    x = (0:(Nx-1))-var_x(j);x = fftshift(x);x(1:(find(x==0)-1))=x(1:(find(x==0)-1))-x(find(x==0)-1)-x(find(x==0)+1);x=ifftshift(x);  
    y = x;
    [X,Y] = meshgrid(x,y);    
    len_k = length(kvec);
    % This simulates an A scan with a 50/50 beam spliter
    E_sc = zeros(Nx,Nx,2048);%scattered field
    
    Er = zeros(Nx,Nx,2048);%reference field
    
    
    %Oringinal E0
    %E0=?
    E0=1;
    % E Scattering spherical wave  
    % Proportion at different scattering unit?
    len_zsc = length(zs);
    temp = zeros(Nx,Nx,len_zsc);
    for ik=1:numel(kvec) 
    % Note, E0 is changable   
        for iz = 1: numel(zs)
            temp(:,:,iz) =exp(-(var_x(j).^2 )/W_0^2) *exp(-1i* kvec(ik) * zs(iz)) * E0./sqrt((X-var_x(j)).^2+Y.^2+zs(iz)^2).*exp(-1i.* kvec(ik).* sqrt((X-var_x(j)).^2+Y.^2+zs(iz).^2));
            %temp(:,:,iz) =ComplexAMP(W_0,2*pi/kvec(ik),[var_x,0,zs(iz)]) * E0./sqrt((X-var_x(j)).^2+Y.^2+zs(iz)^2).*exp(-1i.* kvec(ik).* sqrt((X-var_x(j)).^2+Y.^2+zs(iz).^2));
            
            E_sc(:,:,ik)=E_sc(:,:,ik)+temp(:,:,iz);
        end
    %E_sc(:,:,ik)= fftshift(E_sc(:,:,ik));
    
    %Gaussian Eref
    % E_ref(:,:,ik)=fftshift(100*exp(-(((X).^2 + (Y).^2)/W_0^2)) *exp(-1i* kvec(ik) * zr));
    E_ref(:,:,ik)=(2*exp(-(((X).^2 + (Y).^2)/W_0^2)) *exp(-1i* kvec(ik) * zr));
    
        
    end
    
    fbm = (exp(-(((X).^2 + (Y).^2)/W_0^2)) *exp(-1i* kvec(ik) * 0));
    E_sum = E_sc+E_ref;
    for n = 1:numel(kvec)
        Coupling_Coefficient(n) = sum(sum(E_sum(:,:,n)))*dx*dx;
        Coupling_Coefficient_sc(n) = sum(sum(fbm.*E_sc(:,:,n)))*dx*dx;%Coupling needs to take into account the mode of the fiber
        Coupling_Coefficient_ref(n) = sum(sum(fbm.*E_ref(:,:,n)))*dx*dx;%Coupling needs to take into account the mode of the fiber
    end
    Ik_0 = 1/2*abs(Coupling_Coefficient.^2); %spectral interferogram
    Ik = 2*real(Coupling_Coefficient_sc.*conj(Coupling_Coefficient_ref)); %spectral interferogram
    A = fft(Ik.*tukeywin(numel(Ik)).');
    zbar = (0:(numel(A)-1))/numel(A)/diff(kvec(900:901)/2/pi)/2;
    
    figure(1);clf;
    plot(zbar,abs(A));
    B(j,:)= abs(A);
    j

end

figure(1);clf;
imagesc(zbar,var_x,B)
axis equal
colorbar






