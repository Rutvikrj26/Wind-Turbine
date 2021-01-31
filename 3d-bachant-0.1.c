#include "udf.h"
 
real cl_val [] = {0.   , -0.116, -0.246, -0.375, -0.504, -0.628, -0.756, -0.883,
       -1.007, -1.13 , -1.252, -1.369, -1.484, -1.594, -1.701, -1.804,
       -1.901, -1.992, -2.078, -2.155, -2.222, -2.277, -2.329, -2.37 ,
       -2.402, -2.426, -2.441, -2.445, -2.445, -2.438, -2.424, -2.405,
       -2.382, -2.352, -2.32 , -2.284, -2.245, -2.201, -2.159, -2.116,
       -2.073, -2.028, -1.983, -1.939, -1.895, -1.851, -1.808, -0.798,
       -0.801, -0.803, -0.804, -0.803, -0.798, -0.793, -0.786, -0.78 ,
       -0.773, -0.766, -0.759, -0.752, -0.745, -0.738, -0.731, -0.723,
       -0.717, -0.71 , -0.704, -0.698, -0.693, -0.687, -0.682, -0.678,
       -0.673, -0.669, -0.665, -0.662, -0.659, -0.671, -0.675, -0.679,
       -0.684, -0.688, -0.693, -0.697, -0.701, -0.704, -0.709, -0.713,
       -0.717, -0.721, -0.726, -0.454, -0.46 , -0.465, -0.471, -0.477,
       -0.484, -0.491, -0.498, -0.505, -0.513, -0.521, -0.53 , -0.539,
       -0.548, -0.558, -0.568, -0.579, -0.591, -0.603, -0.615, -0.628,
       -0.642, -0.657, -0.672, -0.688, -0.704, -0.722, -0.74 , -0.759,
       -0.78 , -0.801, -0.823, -0.846, -0.87 , -0.895, -0.921, -0.948,
       -0.977, -1.006, -1.032, -1.059, -1.087, -1.115, -1.144, -1.173,
       -1.203, -1.234, -1.264, -1.295, -1.326, -1.356, -1.386, -1.415,
       -1.444, -1.471, -1.496, -1.518, -1.538, -1.554, -1.568, -1.579,
       -1.584, -1.585, -1.582, -1.575, -1.561, -1.548, -1.525, -1.505,
       -1.48 , -1.452, -1.422, -1.392, -1.36 , -1.328, -1.291, -1.25 ,
       -1.202, -1.146, -1.079, -1.   , -0.911, -0.812, -0.706, -0.593,
       -0.476, -0.356, -0.234, -0.11 , -0.   ,  0.118,  0.243,  0.365,
        0.485,  0.602,  0.715,  0.821,  0.92 ,  1.01 ,  1.089,  1.156,
        1.211,  1.259,  1.3  ,  1.337,  1.369,  1.401,  1.431,  1.461,
        1.488,  1.513,  1.534,  1.556,  1.569,  1.583,  1.59 ,  1.592,
        1.591,  1.586,  1.575,  1.56 ,  1.544,  1.524,  1.501,  1.476,
        1.449,  1.42 ,  1.39 ,  1.36 ,  1.33 ,  1.299,  1.268,  1.238,
        1.207,  1.177,  1.147,  1.118,  1.09 ,  1.062,  1.035,  1.009,
        0.979,  0.951,  0.923,  0.897,  0.872,  0.848,  0.825,  0.803,
        0.782,  0.761,  0.742,  0.724,  0.706,  0.689,  0.673,  0.658,
        0.644,  0.63 ,  0.617,  0.604,  0.592,  0.581,  0.57 ,  0.559,
        0.549,  0.54 ,  0.531,  0.522,  0.514,  0.506,  0.499,  0.492,
        0.485,  0.479,  0.472,  0.466,  0.461,  0.455,  0.45 ,  0.723,
        0.719,  0.714,  0.71 ,  0.705,  0.702,  0.698,  0.694,  0.69 ,
        0.685,  0.681,  0.677,  0.673,  0.66 ,  0.663,  0.667,  0.67 ,
        0.675,  0.679,  0.684,  0.689,  0.694,  0.7  ,  0.706,  0.712,
        0.718,  0.725,  0.732,  0.739,  0.747,  0.754,  0.761,  0.768,
        0.775,  0.782,  0.789,  0.795,  0.801,  0.806,  0.806,  0.806,
        0.804,  0.801,  1.809,  1.851,  1.895,  1.939,  1.983,  2.028,
        2.073,  2.117,  2.161,  2.204,  2.245,  2.284,  2.32 ,  2.353,
        2.382,  2.407,  2.426,  2.44 ,  2.447,  2.448,  2.441,  2.426,
        2.402,  2.37 ,  2.329,  2.28 ,  2.222,  2.155,  2.078,  1.992,
        1.901,  1.804,  1.702,  1.595,  1.485,  1.371,  1.254,  1.134,
        1.012,  0.888,  0.763,  0.637,  0.511,  0.383,  0.256,  0.128,
        0.   };

real cd_val [] = {0.00924, 0.01174, 0.01483, 0.01545, 0.01623, 0.01718, 0.01838,
       0.01994, 0.02178, 0.02386, 0.02619, 0.02893, 0.03198, 0.03533,
       0.0391 , 0.04342, 0.04845, 0.0538 , 0.05997, 0.06629, 0.0735 ,
       0.081  , 0.08978, 0.09916, 0.10956, 0.12025, 0.13017, 0.14374,
       0.15868, 0.17436, 0.19211, 0.20798, 0.227  , 0.24573, 0.26831,
       0.28816, 0.30772, 0.33246, 0.35626, 0.3846 , 0.41256, 0.43707,
       0.46168, 0.48295, 0.51089, 0.55243, 0.56416, 0.63149, 0.67064,
       0.69878, 0.73315, 0.76711, 0.78099, 0.80801, 0.86101, 0.87925,
       0.95641, 0.96622, 0.97586, 1.06363, 1.04518, 1.13003, 1.12784,
       1.12694, 1.21754, 1.25657, 1.26479, 1.26781, 1.32296, 1.33613,
       1.46751, 1.48132, 1.49559, 1.50992, 1.55793, 1.56756, 1.57717,
       1.61246, 1.61534, 1.64907, 1.64719, 1.67664, 1.67039, 1.69582,
       1.72191, 1.70696, 1.72858, 1.75093, 1.72711, 1.74531, 1.76432,
       1.90889, 1.89002, 1.91288, 1.88949, 1.86662, 1.88034, 1.85286,
       1.82581, 1.83019, 1.79872, 1.7984 , 1.7626 , 1.75749, 1.71745,
       1.70753, 1.69755, 1.64915, 1.63439, 1.61966, 1.60535, 1.47358,
       1.45983, 1.40421, 1.40053, 1.39172, 1.35208, 1.2607 , 1.26073,
       1.26195, 1.17623, 1.19376, 1.10505, 1.09432, 1.08337, 1.00529,
       0.98586, 0.9316 , 0.90315, 0.88796, 0.85253, 0.81668, 0.78716,
       0.74648, 0.70669, 0.69269, 0.64886, 0.61859, 0.59509, 0.56792,
       0.54111, 0.51065, 0.47969, 0.45311, 0.4255 , 0.40308, 0.38028,
       0.35487, 0.33322, 0.3114 , 0.29225, 0.27108, 0.25206, 0.23399,
       0.21654, 0.20262, 0.18808, 0.17263, 0.15886, 0.14444, 0.13131,
       0.11828, 0.10593, 0.09347, 0.08208, 0.07132, 0.062  , 0.05403,
       0.04731, 0.04168, 0.03721, 0.03378, 0.03103, 0.02882, 0.02703,
       0.02572, 0.02472, 0.0239 , 0.02326, 0.02017, 0.01767, 0.02016,
       0.02325, 0.02388, 0.02469, 0.02569, 0.02698, 0.02875, 0.03095,
       0.03367, 0.03707, 0.04151, 0.04709, 0.05377, 0.06169, 0.07096,
       0.08168, 0.09303, 0.10546, 0.11777, 0.13078, 0.14388, 0.15829,
       0.17204, 0.18748, 0.20201, 0.21592, 0.23338, 0.25145, 0.27047,
       0.29163, 0.31079, 0.33261, 0.35427, 0.37969, 0.40248, 0.42491,
       0.45252, 0.47911, 0.51008, 0.54055, 0.56736, 0.59455, 0.61806,
       0.64833, 0.69217, 0.70618, 0.74597, 0.78666, 0.8162 , 0.85206,
       0.88749, 0.9027 , 0.93115, 0.98543, 1.00487, 1.08296, 1.09392,
       1.10467, 1.19339, 1.17587, 1.2616 , 1.26039, 1.26037, 1.35177,
       1.39142, 1.40025, 1.40393, 1.45957, 1.47333, 1.60512, 1.61944,
       1.63418, 1.64895, 1.69737, 1.70736, 1.71729, 1.75734, 1.76247,
       1.79829, 1.79861, 1.8301 , 1.82574, 1.8528 , 1.88029, 1.86659,
       1.88947, 1.91287, 1.89003, 1.90892, 1.92848, 1.74537, 1.72717,
       1.751  , 1.72866, 1.70705, 1.72201, 1.69592, 1.6705 , 1.67675,
       1.64732, 1.64921, 1.61549, 1.61262, 1.57734, 1.56774, 1.55811,
       1.5101 , 1.49578, 1.48152, 1.46772, 1.33633, 1.32318, 1.26803,
       1.26501, 1.2568 , 1.21778, 1.12719, 1.12809, 1.13029, 1.04545,
       1.06389, 0.97613, 0.9665 , 0.95669, 0.87954, 0.86131, 0.80831,
       0.78129, 0.76741, 0.73346, 0.6991 , 0.67096, 0.63181, 0.56416,
       0.55243, 0.51089, 0.48295, 0.46168, 0.43707, 0.41256, 0.3846 ,
       0.35626, 0.33246, 0.30772, 0.28816, 0.26831, 0.24573, 0.227  ,
       0.20798, 0.19211, 0.17436, 0.15868, 0.14374, 0.13017, 0.12025,
       0.10956, 0.09916, 0.08978, 0.081  , 0.0735 , 0.06629, 0.05997,
       0.0538 , 0.04845, 0.04342, 0.0391 , 0.03533, 0.03198, 0.02893,
       0.02619, 0.02386, 0.02178, 0.01994, 0.01838, 0.01718, 0.01623,
       0.01545, 0.01483, 0.01174, 0.00925};

real radius = 0.5;
real shaftd = 0.05;
real vel_inf = 10;
real lambda_inf = 1.9;
real rot_vel;
//real timestep = 0.01;
real rho = 1.1225;
real chord_len = 0.14;
real thickness = 0.20*0.14;
real height = 1;
real wings = 3;
real grx = 10;
real gry = 36/3.7;
real chord_approx = 0.14;
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
#define Cd(alpha) cd_val[alpha+180] //replace alpha by alpha_md

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

//    real act_z = 0;
//    real cellid[3];

    rot_vel = w(lambda_inf);

    for (i = 0; i < wings; i++)
    {
        act_theta[i] = theta_inst(rot_vel,time) + i*(2*3.14/wings);
        last_theta[i] = theta_inst(rot_vel,last_time) + i*(2*3.14/wings);
        act_x[i] = x_new(radius, act_theta[i]);
        act_y[i] = y_new(radius, act_theta[i]);
        act_x_mid2[i] = x_new_mid2(radius, act_theta[i]);
        act_y_mid2[i] = y_new_mid2(radius, act_theta[i]);
        act_x_tip[i] = x_new_tip(radius, act_theta[i]);
        act_y_tip[i] = y_new_tip(radius, act_theta[i]);
        act_x_mid[i] = x_new_mid(radius, act_theta[i]);
        act_y_mid[i] = y_new_mid(radius, act_theta[i]);
        act_x_end[i] = x_new_end(radius, act_theta[i]);
        act_y_end[i] = y_new_end(radius, act_theta[i]);

    }

    C_CENTROID(x,c,t);

    if (compshaft(x[0], x[1]))
    {
        source = -0.5*rho*u*u*shaftd*grx*gry/4;
	// printf("Center force value added : %g\n", source );
    }

    for (j = 0; j < wings ; j++)
    {   
        if (fabs(x[2]) <= height/2) 
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
        if (fabs(x[2]) <= height/2) 
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