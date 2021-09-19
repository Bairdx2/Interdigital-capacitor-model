classdef IDC_model < handle
    
    %This is a class definition file that acts as a template
    %for modeling an interdigital capacitor of a certain geometry
    %and with up to one dielectric coating.
    
    properties
        
        L;                          %[m] Length of the IDC fingers
        W;                          %[m] Width of the IDC fingers
        G;                          %[m] Spacing between IDC fingers
        n;                          %[unitless] Metallization ratio
        h;                          %[m] Height of the dielectric layer
        lambda;                     %[m] Twice the length of the sum between the finger width and spacing
        r;                          %[unitless] ration between the dielectric layer height and lambda
        q;
        v2;                         %Second Jacobi function
        v3;                         %Third Jacobi function
        k;                          %Modulus
        z2_ci;                      %Component for measuring interior capacitance
        t2_ci;                      %Component for measuring interior capacitance
        t4_ci;                      %Component for measuring interior capacitance
        k_I;                        %Interior modulus component
        k_Ip;                       %Complimentary interior modulus component
        t3_ce;                      %Component for measuring exterior capacitance
        t4_ce;                      %Component for measuring exterior capacitance
        k_E;                        %Exterior modulus component
        k_Ep;                       %Complimentary exterior modulus component
        k_Ioo;                      %Infinite layer approximation interior modulus
        k_Ioop;                     %Complimentary infinitie layer approximation interior modulus
        k_Eoo;                      %Infinite layer approximation exterior modulus
        k_Eoop;                     %Complimentary infinitie layer approximation exterior modulus
        C_int;                      %Interior capacitance component
        C_ext;                      %Exterior capacitance component
        C_IDC;                      %Total interdigital capacitance
        N;                          %Number of interdigital fingers
        
        
    end
    
    properties (Constant)
        
        eo = 8.854187176E-12;       %[F/m] Permittivity of free space
        er_Al2O3 = 9.43;            %[unitless] Relative permittivity of Al2O3
        er_HfO2 = 19.0;             %[unitless] Relative permittivity of hafnium oxide
        t3_ci = 1;                  %Component for measuring interior capacitance

    end
    
    methods
        
        %The main constructor of this class used to initialize
        %some values when an object is created
        function IDC = IDC_model(Length, Width, Spacing, Height)
            
            IDC.L = Length;
            IDC.W = Width;
            IDC.G = Spacing;
            IDC.h = Height;
            
            %Call to the iteration function which calls all
            %of the other essential functions
            generate_plots(IDC);

            
            
        end
        
        function IDC = plot1(IDC)
            
            %Multiple calls to helper functions that calculate
            %basic values that presumably do not change much
            param_calc(IDC);
            Jacobi2(IDC);
            Jacobi3(IDC);
            k_calc(IDC);
            interior_component_calc(IDC);
            exterior_component_calc(IDC);
            infinite_approx_component_calc(IDC);
            int_ext_capacitance(IDC);
            Total_capacitance1(IDC);
            Cap_vs_nFingers(IDC);

 
        end
        function IDC = plot2(IDC)
           
            Cap_vs_metallization(IDC);
            
        end
        
        function IDC = generate_plots(IDC)
            
            plot1(IDC);
            plot2(IDC);
            %Cap_vs_metallization(IDC);
            
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%Helper functions%%%%%%%%%%%%%%%%%%%%%%%%

function IDC = param_calc(IDC)

    IDC.n = (IDC.W)/(IDC.W + IDC.G);
    IDC.lambda = 2*(IDC.W + IDC.G);
    IDC.r = IDC.h/IDC.lambda;
    IDC.q = exp(-4*pi*IDC.r);

end

function IDC = Jacobi2(IDC)

    %Calculates the second Jacobi function
    v2_sum = 0;
    for w = 0:1000
   
        curr_sum2 = (IDC.q^(w*(w+1)))*cos((2*w+1)*0);
        v2_sum = v2_sum + curr_sum2;
    
    end
    IDC.v2 = 2*(IDC.q^(1/4))*v2_sum;

end

function IDC = Jacobi3(IDC)

    %Calculates the second Jacobi function
    v3_sum = 0;
    for r = 0:1000
    
        curr_sum3 = (IDC.q^((r)^2))*cos(2*r*0);
        v3_sum = v3_sum + curr_sum3;
    
    end
    IDC.v3 = 1 + 2*v3_sum;
    
end

function IDC = k_calc(IDC)

    %Calculates the modulus
    IDC.k = (IDC.v2/IDC.v3)^2;

end

function IDC = interior_component_calc(IDC)

    IDC.z2_ci = ellipticK(IDC.k)*IDC.n;
    [IDC.t2_ci, CN, DN] = ellipj(IDC.z2_ci, IDC.k);
    IDC.t4_ci = 1/IDC.k;
    IDC.k_I = IDC.t2_ci*sqrt(((IDC.t4_ci^2)-1)/((IDC.t4_ci^2)-(IDC.t2_ci^2)));
    IDC.k_Ip = sqrt(1-(IDC.k_I)^2);

end

function IDC = exterior_component_calc(IDC)

    IDC.t3_ce = cosh((pi*(1-IDC.n))/(8*IDC.r));
    IDC.t4_ce = cosh((pi*(1+IDC.n))/(8*IDC.r));
    IDC.k_E = (1/IDC.t3_ce)*sqrt(((IDC.t4_ce^2)-(IDC.t3_ce^2))/((IDC.t4_ce^2)-1));
    IDC.k_Ep = sqrt(1-(IDC.k_E)^2);
    %IDC.k_E = 0;
    %IDC.k_Ep = 0;

end

function IDC = infinite_approx_component_calc(IDC)

    IDC.k_Ioo = sin((pi/2)*IDC.n);
    IDC.k_Ioop = sqrt(1-(IDC.k_Ioo)^2);
    IDC.k_Eoo = (2*sqrt(IDC.n))/(1+IDC.n);
    IDC.k_Eoop = sqrt(1-(IDC.k_Eoo)^2);

end

function IDC = int_ext_capacitance(IDC)

    R1 = ellipticK(IDC.k_Ioo)/ellipticK(IDC.k_Ioop);
    R2 = ellipticK(IDC.k_Eoo)/ellipticK(IDC.k_Eoop);
    R3 = ellipticK(IDC.k_I)/ellipticK(IDC.k_Ip);
    R4 = ellipticK(IDC.k_E)/ellipticK(IDC.k_Ep);
    IDC.C_int = IDC.eo*IDC.L*(R1 + (IDC.er_HfO2-1)*R3 +IDC.er_Al2O3*R1);
    IDC.C_ext = IDC.eo*IDC.L*(R2 + (IDC.er_HfO2-1)*R4 +IDC.er_Al2O3*R2);

end

function IDC = Total_capacitance1(IDC)

    IDC.N = linspace(3,80,100);
    IDC.C_IDC = (IDC.N - 3)*(IDC.C_int/2) + 2*((IDC.C_int * IDC.C_ext)/(IDC.C_int + IDC.C_ext));

end

function IDC = Cap_vs_nFingers(IDC)

    figure;
    plot(IDC.N, IDC.C_IDC*(1E12));
    xlabel('Number of fingers');
    ylabel('Capacitance (pF)');
    grid;

end

function IDC = Cap_vs_metallization(IDC)

    %Calculates the a component of the PR EOS
    N = linspace(20, 100, 5);
    n = linspace(0.1,0.9,100);
    figure();
    xlabel('Metallization ratio (\eta)');
    ylabel('Capacitance (pF)');
    grid;
    hold on;
    for z = 1:5

        N_curr = N(z);
    
        for a = 1:100
        
            Jacobi2(IDC);
            Jacobi3(IDC);
            k_calc(IDC);
            IDC.z2_ci = ellipticK(IDC.k)*n(a);
            [IDC.t2_ci, CN, DN] = ellipj(IDC.z2_ci, IDC.k);
            IDC.t4_ci = 1/IDC.k;
            IDC.k_I = IDC.t2_ci*sqrt(((IDC.t4_ci^2)-1)/((IDC.t4_ci^2)-(IDC.t2_ci^2)));
            IDC.k_Ip = sqrt(1-(IDC.k_I)^2);
            IDC.t3_ce = cosh((pi*(1-n(a)))/(8*IDC.r));
            IDC.t4_ce = cosh((pi*(1+n(a)))/(8*IDC.r));
            IDC.k_E = (1/IDC.t3_ce)*sqrt(((IDC.t4_ce^2)-(IDC.t3_ce^2))/((IDC.t4_ce^2)-1));
            IDC.k_Ep = sqrt(1-(IDC.k_E)^2);
            IDC.k_Ioo = sin((pi/2)*n(a));
            IDC.k_Ioop = sqrt(1-(IDC.k_Ioo)^2);
            IDC.k_Eoo = (2*sqrt(n(a)))/(1+n(a));
            IDC.k_Eoop = sqrt(1-(IDC.k_Eoo)^2);
            R1 = ellipticK(IDC.k_Ioo)/ellipticK(IDC.k_Ioop);
            R2 = ellipticK(IDC.k_Eoo)/ellipticK(IDC.k_Eoop);
            R3 = ellipticK(IDC.k_I)/ellipticK(IDC.k_Ip);
            R4 = ellipticK(IDC.k_E)/ellipticK(IDC.k_Ep);
            IDC.C_int = IDC.eo*IDC.L*(R1 + (IDC.er_HfO2-1)*R3 +IDC.er_Al2O3*R1);
            IDC.C_ext = IDC.eo*IDC.L*(R2 + (IDC.er_HfO2-1)*R4 +IDC.er_Al2O3*R2);
            %C_int = eo*L*(R1 + es*R1);
            %C_ext = eo*L*(R2 + es*R2);
            IDC.C_IDC(a) = (N_curr - 3)*(IDC.C_int/2) + 2*((IDC.C_int * IDC.C_ext)/(IDC.C_int + IDC.C_ext));
        
        end
    plot(n, IDC.C_IDC*(1E12),'color',rand(1,3),'DisplayName',strcat(num2str(N_curr),' fingers'))
end
%plot(n, C_IDC*(1E12));
legend('show')    

end