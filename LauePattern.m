function LauePattern(a,c,lamda_min,dist,phi_o)
%@Nitzan Keynan
%this functions inputs are the folowing:
% dist= distance from the detector
%a,c - information about the lattice, Usually in 10^-10 meter terms
% phi_o=phi_orientation in celcuis, between c&e1
% a is parallel to e2.*e1
 %dist= distance from the detector
%phi_o -celcuis, between c&e1

b=a;
alpha=90 ;%celcuis %alpha=betha=gama
alpha=deg2rad(90);
B_max = 2/lamda_min; %units= meter^-1
phi_o=deg2rad(phi_o);

if lamda_min < 0
    fprintf('not possible')
end

%orientation mat U[A],U[B]
u_a= [a*cos(phi_o+pi/2) a*sin(phi_o+pi/2) 0; 0 0 b; c*cos(phi_o) c*sin(phi_o) 0]';
u_b=(inv(u_a))';

%reciprocal latice parameters
G_reciprocal= ([u_b]')*u_b;
a_recip=sqrt(G_reciprocal(1,1));
b_recip=sqrt(G_reciprocal(2,2));
c_recip=sqrt(G_reciprocal(3,3));

%lattice mat

S= [1 0 0];

% find max val of miller indices- normal_plane B = [h,k,l]
h_max = floor(abs(B_max)/abs(a_recip));
k_max = floor(abs(B_max)/abs(b_recip));
l_max = floor(abs(B_max)/abs(c_recip));

% Show center of radiation (crystal) and points on screen for hkl
point_on_screen_front=[];
point_on_screen_back=[];
for h= -h_max:h_max
    for k=-k_max:k_max
        for l=-l_max:l_max  
            vector_B_hkl = u_b*([h k l]'); % [Bx By Bz]'
            B_hkl= sqrt(dot(vector_B_hkl,vector_B_hkl));
            lamda= (2*(dot(vector_B_hkl,S)))/(B_hkl^2) ;
            r= 1/lamda; %=k0=k1
            if lamda> lamda_min;
                vector_d= lamda*vector_B_hkl - S';
                d= sqrt(dot(vector_d,vector_d));
                f= (-dist)/vector_d(1);
                if vector_d(1)<0 %forward
                    det_x= -(f)*vector_d(3);
                    det_y= (f)*vector_d(2);
                    if det_x < 1*50 && det_y < 1*50 && det_x > - 1*50 && det_y > - 1*50
                        point_on_screen_front=[point_on_screen_front; det_x det_y];
                    end
                elseif vector_d(1)>0  %back
                    det_x= (f)*vector_d(3);
                    det_y= (f)*vector_d(2);
                    if det_x < 1*50 && det_y < 1*50 && det_x > - 1*50 && det_y > - 1*50
                        point_on_screen_back=[point_on_screen_back; det_x det_y];
                    end
                end
            end
        end
    end
end
    subplot(1,2,1);
    plot(point_on_screen_front(:,1),point_on_screen_front(:,2),'red.');
    title('front');
    grid on;
    axis equal;
    subplot(1,2,2);
    plot(point_on_screen_back(:,1),point_on_screen_back(:,2),'red.');
    title('back');
    grid on;
    axis equal;
end
