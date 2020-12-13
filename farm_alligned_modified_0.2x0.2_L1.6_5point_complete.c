#include "udf.h"
 
real cl_val [] = {-0.   , -0.103, -0.232, -0.355, -0.473, -0.6  , -0.721, -0.841,
       -0.953, -1.072, -1.182, -1.286, -1.385, -1.472, -1.556, -1.635,
       -1.696, -1.755, -1.796, -1.83 , -1.852, -1.865, -1.867, -1.86 ,
       -1.838, -1.8  , -1.755, -1.704, -1.649, -1.592, -1.534, -1.471,
       -1.411, -1.351, -1.292, -1.235, -1.18 , -1.127, -1.078, -1.029,
       -0.984, -0.94 , -0.899, -0.861, -0.824, -0.79 , -0.323, -0.32 ,
       -0.317, -0.313, -0.309, -0.305, -0.3  , -0.296, -0.292, -0.287,
       -0.283, -0.279, -0.275, -0.271, -0.267, -0.263, -0.26 , -0.256,
       -0.253, -0.25 , -0.247, -0.244, -0.242, -0.239, -0.237, -0.235,
       -0.233, -0.231, -0.23 , -0.229, -0.229, -0.228, -0.228, -0.228,
       -0.228, -0.229, -0.229, -0.23 , -0.231, -0.232, -0.233, -0.233,
       -0.234, -0.236, -0.148, -0.149, -0.151, -0.153, -0.156, -0.158,
       -0.16 , -0.163, -0.165, -0.168, -0.171, -0.174, -0.177, -0.18 ,
       -0.184, -0.188, -0.192, -0.196, -0.2  , -0.205, -0.21 , -0.215,
       -0.22 , -0.226, -0.232, -0.239, -0.245, -0.253, -0.26 , -0.268,
       -0.277, -0.286, -0.296, -0.306, -0.317, -0.329, -0.341, -0.354,
       -0.368, -0.383, -0.398, -0.415, -0.433, -0.452, -0.472, -0.494,
       -0.517, -0.541, -0.566, -0.594, -0.622, -0.653, -0.685, -0.718,
       -0.753, -0.79 , -0.827, -0.866, -0.907, -0.946, -0.987, -1.026,
       -1.065, -1.101, -1.135, -1.164, -1.189, -1.204, -1.207, -1.206,
       -1.197, -1.181, -1.156, -1.123, -1.082, -1.035, -0.978, -0.912,
       -0.852, -0.883, -0.954, -0.938, -0.877, -0.793, -0.693, -0.585,
       -0.472, -0.354, -0.242, -0.118,  0.   ,  0.123,  0.247,  0.361,
        0.478,  0.592,  0.7  ,  0.8  ,  0.884,  0.945,  0.96 ,  0.89 ,
        0.858,  0.918,  0.985,  1.041,  1.089,  1.129,  1.162,  1.186,
        1.202,  1.211,  1.212,  1.209,  1.193,  1.168,  1.139,  1.105,
        1.069,  1.029,  0.99 ,  0.949,  0.909,  0.869,  0.83 ,  0.792,
        0.755,  0.72 ,  0.686,  0.654,  0.624,  0.595,  0.568,  0.542,
        0.518,  0.495,  0.473,  0.453,  0.434,  0.416,  0.399,  0.384,
        0.369,  0.355,  0.342,  0.329,  0.318,  0.307,  0.296,  0.287,
        0.278,  0.269,  0.261,  0.253,  0.246,  0.239,  0.233,  0.227,
        0.221,  0.215,  0.21 ,  0.205,  0.201,  0.196,  0.192,  0.188,
        0.184,  0.181,  0.177,  0.174,  0.171,  0.168,  0.166,  0.163,
        0.16 ,  0.158,  0.156,  0.154,  0.152,  0.15 ,  0.237,  0.236,
        0.235,  0.234,  0.233,  0.232,  0.231,  0.23 ,  0.23 ,  0.229,
        0.229,  0.228,  0.228,  0.229,  0.229,  0.229,  0.23 ,  0.231,
        0.233,  0.235,  0.237,  0.24 ,  0.242,  0.245,  0.248,  0.25 ,
        0.253,  0.257,  0.26 ,  0.264,  0.267,  0.271,  0.275,  0.279,
        0.284,  0.288,  0.292,  0.297,  0.301,  0.305,  0.31 ,  0.314,
        0.318,  0.321,  0.324,  0.792,  0.827,  0.863,  0.902,  0.943,
        0.986,  1.032,  1.081,  1.132,  1.185,  1.24 ,  1.297,  1.356,
        1.417,  1.478,  1.539,  1.599,  1.658,  1.714,  1.765,  1.811,
        1.849,  1.872,  1.88 ,  1.878,  1.866,  1.844,  1.811,  1.768,
        1.713,  1.648,  1.574,  1.49 ,  1.399,  1.3  ,  1.195,  1.086,
        0.972,  0.854,  0.735,  0.614,  0.492,  0.37 ,  0.247,  0.123,
       -0.   };

real cd_val [] = {0.01709, 0.01782, 0.01864, 0.01566, 0.01997, 0.02125, 0.02298,
       0.02521, 0.02832, 0.03135, 0.03485, 0.03792, 0.04052, 0.04261,
       0.04869, 0.05646, 0.06529, 0.07742, 0.09034, 0.10334, 0.11959,
       0.13804, 0.15233, 0.17413, 0.19587, 0.230071, 0.26041, 0.2759 ,
       0.3171 , 0.37534, 0.40186, 0.44814, 0.48832, 0.51907, 0.5864 ,
       0.62797, 0.67958, 0.78146, 0.79844, 0.86852, 0.94647, 1.00586,
       1.08539, 1.17223, 1.18251, 1.36117, 1.37991, 1.62761, 1.65576,
       1.85314, 1.7991 , 1.99111, 1.95441, 2.0757 , 2.28304, 2.53931,
       2.62769, 2.49196, 2.72737, 3.01291, 2.9749 , 3.23069, 3.06252,
       3.33593, 3.2079 , 3.516  , 3.70913, 3.5234 , 3.82834, 4.18613,
       4.0077 , 3.87963, 4.43597, 4.23024, 4.60947, 4.41683, 4.26098,
       4.67344, 4.52422, 4.40233, 4.88221, 4.76545, 4.67148, 4.59671,
       4.53851, 4.49147, 5.08148, 5.05762, 5.05085, 5.06134, 4.62293,
       5.22123, 5.21028, 5.21644, 5.23993, 4.64954, 4.69544, 4.75263,
       4.82634, 4.91921, 5.0348 , 4.55385, 4.67458, 4.82271, 4.40923,
       4.56389, 4.75551, 4.37574, 4.58098, 4.0242 , 4.1516 , 4.32929,
       3.97084, 3.6652 , 3.85008, 3.65643, 3.34744, 3.4745 , 3.20013,
       3.36734, 3.11054, 3.14765, 2.86102, 2.62448, 2.75896, 2.66951,
       2.41198, 2.20334, 2.08084, 2.1161 , 1.92268, 1.97531, 1.77644,
       1.74664, 1.49751, 1.50724, 1.32635, 1.31366, 1.22438, 1.14246,
       1.08054, 1.00001, 0.92734, 0.90769, 0.80315, 0.74882, 0.7045 ,
       0.63438, 0.60082, 0.55781, 0.50864, 0.47924, 0.41805, 0.3739 ,
       0.35542, 0.32273, 0.28488, 0.26006, 0.23525, 0.21791, 0.19629,
       0.17691, 0.16079, 0.14476, 0.12957, 0.11757, 0.10658, 0.09694,
       0.08615, 0.06459, 0.04729, 0.03964, 0.03529, 0.03159, 0.02914,
       0.02729, 0.02604, 0.02165, 0.02462, 0.02382, 0.02309, 0.02382,
       0.02464, 0.02167, 0.02599, 0.0273 , 0.0291 , 0.03154, 0.03518,
       0.03957, 0.0471 , 0.06426, 0.08567, 0.09644, 0.10607, 0.11701,
       0.12901, 0.1442 , 0.16023, 0.17635, 0.19573, 0.21733, 0.2347 ,
       0.25952, 0.28435, 0.3222 , 0.35488, 0.37338, 0.41752, 0.47872,
       0.50813, 0.55729, 0.60032, 0.63389, 0.70401, 0.74834, 0.80267,
       0.90723, 0.92687, 0.99956, 1.08009, 1.14202, 1.22394, 1.31323,
       1.32593, 1.50683, 1.4971 , 1.74624, 1.77605, 1.97493, 1.9223 ,
       2.11573, 2.08048, 2.20299, 2.41164, 2.66918, 2.75864, 2.62416,
       2.86072, 3.14736, 3.11025, 3.36706, 3.19986, 3.47424, 3.34719,
       3.65618, 3.84985, 3.66498, 3.97063, 4.32909, 4.1514 , 4.02402,
       4.58081, 4.37558, 4.75535, 4.56375, 4.40909, 4.82259, 4.67447,
       4.55375, 5.03471, 4.91913, 4.82627, 4.75257, 4.69539, 4.6495 ,
       5.2399 , 5.21643, 5.21027, 5.22123, 4.46252, 5.06138, 5.05089,
       5.05767, 5.08153, 4.49153, 4.53857, 4.59678, 4.67155, 4.76553,
       4.88229, 4.40243, 4.52431, 4.67355, 4.26109, 4.41695, 4.60958,
       4.23037, 4.4361 , 3.87976, 4.00784, 4.18627, 3.82849, 3.52355,
       3.70929, 3.51615, 3.20806, 3.3361 , 3.06269, 3.23087, 2.97508,
       3.01309, 2.72756, 2.49215, 2.62788, 2.53951, 2.28324, 2.07591,
       1.95462, 1.99132, 1.79932, 1.85336, 1.65598, 1.62784, 1.38014,
       1.36117, 1.18251, 1.17223, 1.08538, 1.00585, 0.94647, 0.86852,
       0.79844, 0.78146, 0.67958, 0.62797, 0.5864 , 0.51907, 0.48831,
       0.44814, 0.40185, 0.37534, 0.3171 , 0.27591, 0.26041, 0.23071,
       0.19586, 0.17413, 0.15233, 0.13808, 0.11959, 0.10335, 0.09034,
       0.07743, 0.06528, 0.05646, 0.04864, 0.04256, 0.04052, 0.03793,
       0.03486, 0.03131, 0.02836, 0.02522, 0.02298, 0.02122, 0.02001,
       0.01564, 0.01862, 0.01782, 0.01709};

real radius = 1.4;
real shaftd = 0.1;
real vel_inf = 10;
real lambda_inf = 1.6;
real rot_vel;
//real timestep = 0.01;
real rho = 1.1225;
real chord_len = 0.42;
real thickness = 0.15*0.42;
int wings = 3;
real turb_dist_x = 30;
real turb_dist_y = 30;
int turb_x = 3;
int turb_y = 4;
int tot_turb = turb_x*turb_y
int tot_wings = turb_x*turb_y*wings;
real turb_loc_x[] = {15,15,15,15,45,45,45,45,75,75,75,75};
real turb_loc_y[] = {15,45,75,105,15,45,75,105,15,45,75,105};
// domain size = 150*110
real grx = 5;
real gry = 5;
real chord_approx = 0.42;
real md[] = {0.2,0.2,0.2,0.2,0.2};
real sd[] = {0.2,0.1};
/* general use functions defined here */

#define sign(val) val==0 ? 0 : val/fabs(val) /* can create a singularity issue */

/* Analytical Solution of Theta Defined here*/

#define theta_inst(w,t) w*t
#define x_new(r,theta) r*cos(theta)
#define y_new(r,theta) r*sin(theta) 
#define x_new_tip(r,theta) r*cos(theta) - (chord_approx/4)*sin(theta)
#define y_new_tip(r,theta) r*sin(theta) + (chord_approx/4)*cos(theta)
#define x_new_mid(r,theta) r*cos(theta) + (chord_approx/4)*sin(theta)
#define y_new_mid(r,theta) r*sin(theta) - (chord_approx/4)*cos(theta)
#define x_new_mid2(r,theta) r*cos(theta) + (2*chord_approx/4)*sin(theta)
#define y_new_mid2(r,theta) r*sin(theta) - (2*chord_approx/4)*cos(theta)
#define x_new_end(r,theta) r*cos(theta) + (3*chord_approx/4)*sin(theta)
#define y_new_end(r,theta) r*sin(theta) - (3*chord_approx/4)*cos(theta)

#define beta(theta, u, v) (theta + atan(v/u))      /* can create a singularity issue */
#define alpha(lambda, beta) acos((lambda - sin(beta))/sqrt(pow((cos(beta)),2) + pow((lambda - sin(beta)),2)))   
#define v_net(u,v) sqrt(u*u + v*v)
#define v_rel(v_net, lambda, beta) (v_net*pow((pow((lambda-sin(beta)),2)+pow((cos(beta)),2)),0.5))
#define lambda(v_net) (radius*rot_vel/v_net)
#define w(lbd) vel_inf*lbd/radius

#define gamma_l  (1.4 - 6*(0.06 - (thickness/chord_len)))
#define gamma_d  (1 - 2.5*(0.06 - (thickness/chord_len))) 
#define alpha_dot(alpha_n, alpha_l, timestep) ((alpha_n - alpha_l)/timestep)  /* can create a singularity issue */
#define kappa(dot_sign) 0.75 + 0.25*dot_sign
#define del_l(alpha,kappa,alpha_dot, relative_vel,dot_sign) (gamma_l*kappa*sqrt(fabs(0.5*chord_len*alpha_dot/relative_vel))*dot_sign)  /* can create a singularit☻y issue */
#define del_d(alpha,kappa,alpha_dot, relative_vel,dot_sign) (gamma_d*kappa*sqrt(fabs(0.5*chord_len*alpha_dot/relative_vel))*dot_sign)  /* can create a singularity issue */

#define Cl(alpha) cl_val[alpha+180] //replace alpha by alpha_ml
#define Cd(alpha) cd_val[alpha+180] //replacce alpha by alpha_md

#define Ct(C_lm,C_dm,alpha)  (fabs(C_lm)*sin(alpha)-fabs(C_dm)*cos(alpha))
#define Cn(C_lm,C_dm,alpha)  (fabs(C_lm)*cos(alpha)+fabs(C_dm)*sin(alpha))
/* Calculating Forces from Coefficients */

#define Fn(relative_vel, Cn) 0.5*rho*chord_len*(relative_vel*relative_vel)*Cn*grx*gry
#define Ft(relative_vel, Ct) 0.5*rho*chord_len*(relative_vel*relative_vel)*Ct*grx*gry

#define Fx(relative_vel, Cx) 0.5*rho*chord_len*(relative_vel*relative_vel)*Cx*grx*gry
#define Fy(relative_vel, Cy) 0.5*rho*chord_len*(relative_vel*relative_vel)*Cy*grx*gry

/* Comparision Functions */

#define compx(x_cent,x_act) (fabs(x_act - x_cent) < 0.5*(pow(grx,-1))) 
#define compy(y_cent,y_act) (fabs(y_act - y_cent) < 0.5*(pow(gry,-1))) // && (y_act =< y_cent + 0.5*(pow(gry,-1)))) 
#define compshaft(x_cent,y_cent) (fabs(sqrt(x_cent*x_cent + y_cent*y_cent)) <= 0.5*shaftd) 
//#define comprad(x_cent,y_cent) (fabs((sqrt(x_cent*x_cent + x_cent*x_cent)) - radius) <= tol) 
//#define comparg(x_cent, y_cent, req_theta) (fabs((180/3.14)*(atan(y_cent/x_cent) - req_theta)) < 10)
#define compxrel(x_cent,x_act) (fabs(x_act - x_cent) >= 0.5*(pow(grx,-1)) && fabs(x_act - x_cent) < 1.5*(pow(grx,-1)))
#define compyrel(y_cent,y_act) (fabs(y_act - y_cent) >= 0.5*(pow(gry,-1)) && fabs(y_act - y_cent) < 1.5*(pow(gry,-1)))

/* Actual Source Commands start After here */

 DEFINE_SOURCE(xmom_source, c, t, dS, eqn)
 {
    real source = 0;
    real x[ND_ND];
    real time = CURRENT_TIME;
    real last_time = PREVIOUS_TIME;
    real timestep = (time - last_time);
    real x_vel = C_U(c,t);
    real y_vel = C_V(c,t);

    real x_vel_last = C_U_M1(c,t);
    real y_vel_last = C_V_M1(c,t);
     
    real u = x_vel;
    real v = y_vel;

    real vnet_act = v_net(u,v);

    real u_last = x_vel_last;
    real v_last = y_vel_last;
    real v_net_last = v_net(u_last, v_last);

    int thet,i,j;

    real act_theta [tot_wings];
    real last_theta [tot_wings];
    real act_x [tot_wings];
    real act_y [tot_wings];
    real act_x_tip [tot_wings];
    real act_y_tip [tot_wings];
    real act_x_mid [tot_wings];
    real act_y_mid [tot_wings];
    real act_x_mid2 [tot_wings];
    real act_y_mid2 [tot_wings];
    real act_x_end [tot_wings];
    real act_y_end [tot_wings]; 

//    real act_z = 0;
//    real cellid[3];

    rot_vel = w(lambda_inf);

    for (thet = 0; thet <= tot_wings; (thet=thet+3))
    {
        act_theta[thet] = theta_inst(rot_vel,time);
        act_theta[thet+1] = theta_inst(rot_vel,time) + (2*3.14/wings);
        act_theta[thet+2] = theta_inst(rot_vel,time) + 2*(2*3.14/wings);
        last_theta[thet] = theta_inst(rot_vel,last_time);
        last_theta[thet+1] = theta_inst(rot_vel,last_time) + (2*3.14/wings);
        last_theta[thet+2] = theta_inst(rot_vel,last_time) + 2*(2*3.14/wings);
    }

    for (i = 0; i < tot_wings; i++)
    {
        act_x[i] = turb_loc_x[i] + x_new(radius, act_theta[i]);
        act_y[i] = turb_loc_y[i] + y_new(radius, act_theta[i]);
        act_x_mid2[i] = turb_loc_x[i] + x_new_mid2(radius, act_theta[i]);
        act_y_mid2[i] = turb_loc_y[i] + y_new_mid2(radius, act_theta[i]);
        act_x_tip[i] = turb_loc_x[i] + x_new_tip(radius, act_theta[i]);
        act_y_tip[i] = turb_loc_y[i] + y_new_tip(radius, act_theta[i]);
        act_x_mid[i] = turb_loc_x[i] + x_new_mid(radius, act_theta[i]);
        act_y_mid[i] = turb_loc_y[i] + y_new_mid(radius, act_theta[i]);
        act_x_end[i] = turb_loc_x[i] + x_new_end(radius, act_theta[i]);
        act_y_end[i] = turb_loc_y[i] + y_new_end(radius, act_theta[i]);

    }

    C_CENTROID(x,c,t);

    if (compshaft(x[0], x[1]))
    {
        source = -0.5*rho*u*u*shaftd*grx*gry/4;
	// printf("Center force value added : %g\n", source );
    }

    for (j = 0; j < wings ; j++)
    {
        if (compx(x[0], act_x[j]) && compy(x[1], act_y[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);

                FILE *data;
                data=fopen("data.txt","a");
                fprintf(data,"%d,%g,%g,%g,%g,%g,%g,%g,%g\n", j,(180/3.14)*act_theta[j],alpha_act,alpha_dot_act,alpha_ml_act,alpha_md_act,v_rel_act, Ct_act, Ft_act);
                fclose(data);

                if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                    {
                        source = sd[0]*md[1]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j])) + source;
                    }
                else
                    {
                        source = sd[0]*md[1]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                    }
            

                

//         printf("%g\n", source);    
//         printf("Focal Coordinates Matched\n");
//         printf("theta -> %g\n", (180/3.14)*act_theta[j]);
//         printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
//         printf("u: %g , v: %g\n", u, v );
//         printf("u-last: %g , v-last: %g\n", u_last, v_last);
//         printf("rel velocity : %g\n", v_rel_act);
//         printf("lbd_act : %g, beta_act : %g\n", lbd_act, beta_act);
//         printf("alpha : %g, alpha_last : %g, alpha_dot : %g\n", alpha_act, alpha_last, alpha_dot_act);
//         printf("alpha_ml : %g, alpha_md : %g\n", alpha_ml_act, alpha_md_act);
//         printf("kappa : %g, sign : %g\n", kappa_act, sign(alpha_dot_act));
//         printf("del_l : %g, del_d : %g\n", del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign), del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));         printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
//         printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
//         printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
//         printf("force-x : %g\n", -source);
//         printf("\n");    
        dS[eqn] = 0;
        }
        if (compxrel(x[0], act_x[j]) && compyrel(x[1], act_y[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[1]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j])) + source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[1]*md[1]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("Focal related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }
        if (compx(x[0], act_x[j]) && compyrel(x[1], act_y[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[1]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j])) + source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[1]*md[1]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("Focal related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }
        if (compxrel(x[0], act_x[j]) && compy(x[1], act_y[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[1]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j])) + source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[1]*md[1]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("Focal related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }        

        if (compx(x[0], act_x_tip[j]) && compy(x[1], act_y_tip[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[0]*md[0]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[0]*md[0]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("Tip Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }
        if (compxrel(x[0], act_x_tip[j]) && compyrel(x[1], act_y_tip[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[0]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[1]*md[0]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("Tip related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }        
        if (compxrel(x[0], act_x_tip[j]) && compy(x[1], act_y_tip[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[0]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[1]*md[0]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("Tip related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }        
        if (compx(x[0], act_x_tip[j]) && compyrel(x[1], act_y_tip[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[0]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[1]*md[0]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("Tip related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }        

        if (compx(x[0], act_x_mid[j]) && compy(x[1], act_y_mid[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[0]*md[2]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[0]*md[2]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("Mid pt Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }
        if (compxrel(x[0], act_x_mid[j]) && compyrel(x[1], act_y_mid[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[2]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[1]*md[2]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("Mid related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }
        if (compx(x[0], act_x_mid[j]) && compyrel(x[1], act_y_mid[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[2]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[1]*md[2]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("Mid related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }
        if (compxrel(x[0], act_x_mid[j]) && compy(x[1], act_y_mid[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[2]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[1]*md[2]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("Mid related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }

        if (compx(x[0], act_x_mid2[j]) && compy(x[1], act_y_mid2[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[0]*md[2]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[0]*md[2]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("Mid pt Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }
        if (compxrel(x[0], act_x_mid2[j]) && compyrel(x[1], act_y_mid2[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[2]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[1]*md[2]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("Mid related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }
        if (compx(x[0], act_x_mid2[j]) && compyrel(x[1], act_y_mid2[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[2]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[1]*md[2]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("Mid related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }
        if (compxrel(x[0], act_x_mid2[j]) && compy(x[1], act_y_mid2[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[2]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[1]*md[2]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("Mid related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }

        if (compx(x[0], act_x_end[j]) && compy(x[1], act_y_end[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[0]*md[3]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[0]*md[3]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("End pt Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }
        if (compxrel(x[0], act_x_end[j]) && compyrel(x[1], act_y_end[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[3]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[1]*md[3]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("End pt related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }
        if (compx(x[0], act_x_end[j]) && compyrel(x[1], act_y_end[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[3]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[1]*md[3]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("End pt related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }
        if (compxrel(x[0], act_x_end[j]) && compy(x[1], act_y_end[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[3]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[1]*md[3]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("End pt related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }

    }

    return source;
 } 

 DEFINE_SOURCE(ymom_source, c, t, dS, eqn)
 {
    real source = 0;
    real x[ND_ND];
    real time = CURRENT_TIME;
    real last_time = PREVIOUS_TIME;
    real timestep = time - last_time;
    real x_vel = C_U(c,t);
    real y_vel = C_V(c,t);


    real x_vel_last = C_U_M1(c,t);
    real y_vel_last = C_V_M1(c,t);

    real u = x_vel;
    real v = y_vel;

    real vnet_act = v_net(u,v);
     
    real u_last = x_vel_last;
    real v_last = y_vel_last;
    real v_net_last = v_net(u_last, v_last);

    int i,j;

    real act_theta [3];
    real last_theta [3];
    real act_x [3];
    real act_y [3];
    real act_x_tip [3];
    real act_y_tip [3];
    real act_x_mid [3];
    real act_y_mid [3];    
    real act_x_mid2 [3];
    real act_y_mid2 [3];    
    real act_x_end [3];
    real act_y_end [3];

    rot_vel = w(lambda_inf);

    for (i = 0; i < wings; i++)
    {
        act_theta[i] = theta_inst(rot_vel,time) + i*(2*3.14/wings);
        last_theta[i] = theta_inst(rot_vel, last_time) + i*(2*3.14/wings); 
        act_x[i] = x_new(radius, act_theta[i]);
        act_y[i] = y_new(radius, act_theta[i]);
        act_x_tip[i] = x_new_tip(radius, act_theta[i]);
        act_y_tip[i] = y_new_tip(radius, act_theta[i]);
        act_x_mid[i] = x_new_mid(radius, act_theta[i]);
        act_y_mid[i] = y_new_mid(radius, act_theta[i]);
        act_x_mid2[i] = x_new_mid(radius, act_theta[i]);
        act_y_mid2[i] = y_new_mid(radius, act_theta[i]);
        act_x_end[i] = x_new_end(radius, act_theta[i]);
        act_y_end[i] = y_new_end(radius, act_theta[i]);
    }

    C_CENTROID(x,c,t);

    for (j = 0; j < wings ; j++)
    {
        if (compx(x[0], act_x[j]) && compy(x[1], act_y[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);

            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[0]*md[1]*(Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source;
          //  source = Fy(v_rel_act, Cy_act);
                }
            else
                {
                    source = sd[0]*md[1]*(-Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source; 
                }
            
        //   // printf("%g\n", source);    
            dS[eqn] = 0;
        // printf("Focal Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-y : %g\n", -source);
        // printf("\n");    
        }
        if (compxrel(x[0], act_x[j]) && compyrel(x[1], act_y[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);

            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[1]*(Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source;
          //  source = Fy(v_rel_act, Cy_act);
                }
            else
                {
                    source = sd[1]*md[1]*(-Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source; 
                }
            
        //   // printf("%g\n", source);    
            dS[eqn] = 0;
        // printf("focal related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-y : %g\n", -source);
        // printf("\n");    
        }
        if (compx(x[0], act_x[j]) && compyrel(x[1], act_y[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);

            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[1]*(Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source;
          //  source = Fy(v_rel_act, Cy_act);
                }
            else
                {
                    source = sd[1]*md[1]*(-Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source; 
                }
            
        //   // printf("%g\n", source);    
            dS[eqn] = 0;
        // printf("focal related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-y : %g\n", -source);
        // printf("\n");    
        }
        if (compxrel(x[0], act_x[j]) && compy(x[1], act_y[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);

            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[1]*(Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source;
          //  source = Fy(v_rel_act, Cy_act);
                }
            else
                {
                    source = sd[1]*md[1]*(-Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source; 
                }
            
        //   // printf("%g\n", source);    
            dS[eqn] = 0;
        // printf("focal related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-y : %g\n", -source);
        // printf("\n");    
        }

        if (compx(x[0], act_x_tip[j]) && compy(x[1], act_y_tip[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);

            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[0]*md[0]*(Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source;
          //  source = Fy(v_rel_act, Cy_act);
                }
            else
                {
                    source = sd[0]*md[0]*(-Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source; 
                }
            
        //   // printf("%g\n", source);    
            dS[eqn] = 0;
        // printf(" Tip Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-y : %g\n", -source);
        // printf("\n");    
        }
        if (compxrel(x[0], act_x_tip[j]) && compyrel(x[1], act_y_tip[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);

            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[0]*(Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source;
          //  source = Fy(v_rel_act, Cy_act);
                }
            else
                {
                    source = sd[1]*md[0]*(-Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source; 
                }
            
        //   // printf("%g\n", source);    
            dS[eqn] = 0;
        // printf(" Tip related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-y : %g\n", -source);
        // printf("\n");    
        }
        if (compx(x[0], act_x_tip[j]) && compyrel(x[1], act_y_tip[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);

            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[0]*(Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source;
          //  source = Fy(v_rel_act, Cy_act);
                }
            else
                {
                    source = sd[1]*md[0]*(-Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source; 
                }
            
        //   // printf("%g\n", source);    
            dS[eqn] = 0;
        // printf(" Tip related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-y : %g\n", -source);
        // printf("\n");    
        }
        if (compxrel(x[0], act_x_tip[j]) && compy(x[1], act_y_tip[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);

            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[0]*(Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source;
          //  source = Fy(v_rel_act, Cy_act);
                }
            else
                {
                    source = sd[1]*md[0]*(-Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source; 
                }
            
        //   // printf("%g\n", source);    
            dS[eqn] = 0;
        // printf(" Tip related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-y : %g\n", -source);
        // printf("\n");    
        }

        if (compx(x[0], act_x_mid[j]) && compy(x[1], act_y_mid[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);

            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[0]*md[2]*(Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source;
          //  source = Fy(v_rel_act, Cy_act);
                }
            else
                {
                    source = sd[0]*md[2]*(-Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source; 
                }
            
        //   // printf("%g\n", source);    
            dS[eqn] = 0;
        // printf("Mid Pt Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-y : %g\n", -source);
        // printf("\n");    
        }
        if (compxrel(x[0], act_x_mid[j]) && compyrel(x[1], act_y_mid[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);

            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[2]*(Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source;
          //  source = Fy(v_rel_act, Cy_act);
                }
            else
                {
                    source = sd[1]*md[2]*(-Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source; 
                }
            
        //   // printf("%g\n", source);    
            dS[eqn] = 0;
        // printf("Mid Pt related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-y : %g\n", -source);
        // printf("\n");    
        }
        if (compx(x[0], act_x_mid[j]) && compyrel(x[1], act_y_mid[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);

            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[2]*(Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source;
          //  source = Fy(v_rel_act, Cy_act);
                }
            else
                {
                    source = sd[1]*md[2]*(-Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source; 
                }
            
        //   // printf("%g\n", source);    
            dS[eqn] = 0;
        // printf("Mid Pt related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-y : %g\n", -source);
        // printf("\n");    
        }
        if (compxrel(x[0], act_x_mid[j]) && compy(x[1], act_y_mid[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);

            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[2]*(Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source;
          //  source = Fy(v_rel_act, Cy_act);
                }
            else
                {
                    source = sd[1]*md[2]*(-Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source; 
                }
            
        //   // printf("%g\n", source);    
            dS[eqn] = 0;
        // printf("Mid Pt related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-y : %g\n", -source);
        // printf("\n");    
        }

        if (compx(x[0], act_x_mid2[j]) && compy(x[1], act_y_mid2[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[0]*md[2]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[0]*md[2]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("Mid pt Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }
        if (compxrel(x[0], act_x_mid2[j]) && compyrel(x[1], act_y_mid2[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[2]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[1]*md[2]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("Mid related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }
        if (compx(x[0], act_x_mid2[j]) && compyrel(x[1], act_y_mid2[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[2]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[1]*md[2]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("Mid related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }
        if (compxrel(x[0], act_x_mid2[j]) && compy(x[1], act_y_mid2[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);
         
            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[2]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source;
//            source = Fx(v_rel_act, Cx_act);
                }
            else
                {
                    source = sd[1]*md[2]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                }
            
//        // printf("%g\n", source);    
        // printf("Mid related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }

        if (compx(x[0], act_x_end[j]) && compy(x[1], act_y_end[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);

            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[0]*md[3]*(Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source;
          //  source = Fy(v_rel_act, Cy_act);
                }
            else
                {
                    source = sd[0]*md[3]*(-Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source; 
                }
            
        //   // printf("%g\n", source);    
            dS[eqn] = 0;
        // printf("End Pt Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);        
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-y : %g\n", -source);
        // printf("\n");    
        }
        if (compxrel(x[0], act_x_end[j]) && compyrel(x[1], act_y_end[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);

            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[3]*(Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source;
          //  source = Fy(v_rel_act, Cy_act);
                }
            else
                {
                    source = sd[1]*md[3]*(-Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source; 
                }
            
        //   // printf("%g\n", source);    
            dS[eqn] = 0;
        // printf("End pt related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-y : %g\n", -source);
        // printf("\n");    
        }
        if (compx(x[0], act_x_end[j]) && compyrel(x[1], act_y_end[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);

            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[3]*(Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source;
          //  source = Fy(v_rel_act, Cy_act);
                }
            else
                {
                    source = sd[1]*md[3]*(-Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source; 
                }
            
        //   // printf("%g\n", source);    
            dS[eqn] = 0;
        // printf("End pt related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-y : %g\n", -source);
        // printf("\n");    
        }
        if (compxrel(x[0], act_x_end[j]) && compy(x[1], act_y_end[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_last = lambda(v_net_last);            
            real beta_act = beta(act_theta[j], u, v);
            real beta_last = beta(last_theta[j], u_last, v_last);
            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            real v_rel_last = v_rel(v_net_last, lbd_last, beta_last);                            
            real alpha_act = alpha(lbd_act, beta_act);                
            real alpha_last = alpha(lbd_last, beta_last);                
            real alpha_dot_act = alpha_dot(alpha_act, alpha_last,timestep);
            real alpha_dot_sign = sign(alpha_dot_act);
            real kappa_act = kappa(alpha_dot_sign);                 
            real alpha_ml_act = 57.32*(alpha_act - del_l(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real alpha_md_act = 57.32*(alpha_act - del_d(alpha_act, kappa_act, alpha_dot_act, v_rel_act,alpha_dot_sign));
            real Cl_act = (alpha_act*57.32/alpha_ml_act)*Cl((int) alpha_ml_act);
            real Cd_act = Cd((int) alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
                  
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);

            if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                {
                    source = sd[1]*md[3]*(Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source;
          //  source = Fy(v_rel_act, Cy_act);
                }
            else
                {
                    source = sd[1]*md[3]*(-Fn_act * sin (act_theta[j]) - Ft_act * cos (act_theta[j])) + source; 
                }
            
        //   // printf("%g\n", source);    
            dS[eqn] = 0;
        // printf("End pt related Coordinates Matched\n");
        // printf("theta -> %g\n", (180/3.14)*act_theta[j]);
        // printf("x-coordinate: %g , y-coordinate: %g\n", x[0] , x[1]);
        // printf("u: %g , v: %g\n", u, v );
        // printf("u-last: %g , v-last: %g\n", u_last, v_last);
        // printf("net velocity : %g\n", vnet_act);
        // printf("net velocity : %g\n", v_net_last);
        // printf("alpha : %d, alpha_last : %d, alpha_dot : %g\n", alpha_act, alpha_act_last, alpha_dot_act);
        // printf("alpha_ml : %d, alpha_md : %d\n", alpha_ml_act, alpha_md_act);
        // printf("Cl_m : %g, Cd_m : %g\n", Cl_act, Cd_act);
        // printf("Ct : %g, Cn : %g\n", Ct_act, Cn_act);
        // printf("Ft : %g, Fn : %g\n", Ft_act, Fn_act);
        // printf("force-y : %g\n", -source);
        // printf("\n");    
        }

    }

    return source;
 } 

DEFINE_DELTAT(timestepping, d)
{
    real CFL = 0.5;
    real max_vel = 1;
    real time_step;

    Thread *t;
    thread_loop_c(t, d) /*loops over all cell threads in domain*/
    {                        
            cell_t c;
            begin_c_loop(c, t)    /* loops over cells in a cell thread  */
            {               
                real u = C_U(c,t);
                real v = C_V(c,t);          
                real vel = v_net(u,v);
                if (vel > max_vel)
                {
                    max_vel = vel;
                }
            }                         
            end_c_loop(c, t)
    }
    time_step = CFL/(grx*max_vel);
    return time_step;
}