#include "udf.h"
 
real cl_val [] = {0.   ,  0.123,  0.246,  0.369,  0.492,  0.604,  0.718,  0.825,
        0.924,  1.009,  1.073,  1.104,  1.098,  1.071,  1.031,  1.047,
        1.084,  1.122,  1.155,  1.182,  1.201,  1.213,  1.219,  1.219,
        1.212,  1.2  ,  1.184,  1.163,  1.14 ,  1.113,  1.085,  1.055,
        1.02 ,  0.981,  0.942,  0.903,  0.865,  0.828,  0.792,  0.758,
        0.725,  0.694,  0.664,  0.635,  0.608,  0.583,  0.559,  0.536,
        0.514,  0.494,  0.475,  0.457,  0.439,  0.423,  0.408,  0.394,
        0.38 ,  0.367,  0.355,  0.344,  0.333,  0.323,  0.314,  0.305,
        0.296,  0.288,  0.28 ,  0.273,  0.266,  0.26 ,  0.254,  0.248,
        0.243,  0.237,  0.232,  0.228,  0.223,  0.219,  0.215,  0.211,
        0.207,  0.204,  0.201,  0.198,  0.195,  0.192,  0.189,  0.187,
        0.184,  0.182,  0.301,  0.3  ,  0.298,  0.297,  0.295,  0.294,
        0.292,  0.291,  0.289,  0.288,  0.287,  0.285,  0.284,  0.282,
        0.281,  0.28 ,  0.278,  0.28 ,  0.282,  0.284,  0.287,  0.289,
        0.292,  0.295,  0.298,  0.302,  0.305,  0.309,  0.313,  0.317,
        0.322,  0.326,  0.331,  0.335,  0.34 ,  0.345,  0.35 ,  0.355,
        0.36 ,  0.364,  0.369,  0.373,  0.378,  0.381,  0.384,  0.936,
        0.974,  1.015,  1.058,  1.103,  1.151,  1.201,  1.253,  1.307,
        1.363,  1.42 ,  1.479,  1.539,  1.599,  1.652,  1.698,  1.741,
        1.782,  1.819,  1.851,  1.877,  1.897,  1.909,  1.913,  1.907,
        1.891,  1.865,  1.829,  1.782,  1.724,  1.657,  1.58 ,  1.495,
        1.402,  1.302,  1.196,  1.085,  0.97 ,  0.853,  0.734,  0.614,
        0.492,  0.369,  0.246,  0.123, -0.   , -0.123, -0.246, -0.369,
       -0.492, -0.614, -0.736, -0.857, -0.979, -1.102, -1.23 , -1.365,
       -1.515, -1.678, -1.86 , -2.016, -2.146, -2.257, -2.356, -2.437,
       -2.505, -2.561, -2.601, -2.627, -2.642, -2.648, -2.639, -2.625,
       -2.6  , -2.567, -2.529, -2.487, -2.43 , -2.362, -2.293, -2.221,
       -2.151, -2.081, -2.013, -1.947, -1.879, -1.817, -1.758, -1.7  ,
       -1.644, -1.589, -1.539, -1.491, -1.444, -1.401, -1.359, -1.32 ,
       -1.282, -1.246, -1.213, -1.181, -1.15 , -1.122, -1.095, -1.069,
       -1.044, -1.021, -1.   , -0.979, -0.96 , -0.942, -0.925, -0.909,
       -0.895, -0.881, -0.867, -0.855, -0.844, -0.833, -0.823, -0.814,
       -0.806, -0.798, -0.791, -0.784, -0.779, -0.773, -0.769, -0.765,
       -0.761, -0.758, -0.756, -0.754, -0.753, -0.752, -0.515, -0.517,
       -0.52 , -0.523, -0.526, -0.531, -0.536, -0.542, -0.549, -0.556,
       -0.564, -0.573, -0.582, -0.592, -0.603, -0.614, -0.627, -0.635,
       -0.644, -0.653, -0.663, -0.674, -0.685, -0.697, -0.71 , -0.724,
       -0.739, -0.754, -0.771, -0.788, -0.807, -0.827, -0.848, -0.87 ,
       -0.894, -0.919, -0.946, -0.974, -1.004, -1.035, -1.068, -1.104,
       -1.141, -1.18 , -1.222, -0.936, -0.974, -1.015, -1.058, -1.103,
       -1.151, -1.201, -1.253, -1.307, -1.363, -1.42 , -1.479, -1.539,
       -1.599, -1.652, -1.698, -1.741, -1.782, -1.819, -1.851, -1.877,
       -1.897, -1.909, -1.913, -1.907, -1.891, -1.865, -1.829, -1.782,
       -1.724, -1.657, -1.58 , -1.495, -1.402, -1.302, -1.196, -1.085,
       -0.97 , -0.853, -0.734, -0.614, -0.492, -0.369, -0.246, -0.123,
        0.   };

real cd_val [] = {1.29700e-02,  1.30600e-02,  1.33600e-02,  1.38700e-02,
        1.46800e-02,  1.33500e-02,  1.72100e-02,  1.87600e-02,
        2.05200e-02,  2.25600e-02,  2.56400e-02,  3.03900e-02,
        3.87400e-02,  5.15000e-02,  6.99500e-02,  8.34800e-02,
        9.40100e-02,  1.03800e-01,  1.12940e-01,  1.23820e-01,
        1.34320e-01,  1.46670e-01,  1.57460e-01,  1.73350e-01,
        1.86290e-01,  2.00570e-01,  2.12840e-01,  2.29770e-01,
        2.53000e-01,  2.63200e-01,  2.80070e-01,  3.00140e-01,
        3.24430e-01,  3.53120e-01,  3.68080e-01,  3.93340e-01,
        4.26030e-01,  4.68210e-01,  4.67080e-01,  5.16850e-01,
        5.40940e-01,  5.66410e-01,  5.99780e-01,  6.27680e-01,
        6.81910e-01,  6.95170e-01,  7.27820e-01,  8.06760e-01,
        8.03630e-01,  8.28770e-01,  8.62670e-01,  9.25010e-01,
        1.00321e+00,  1.04165e+00,  1.09413e+00,  1.05220e+00,
        1.06482e+00,  1.13518e+00,  1.18738e+00,  1.14285e+00,
        1.25130e+00,  1.33093e+00,  1.25839e+00,  1.41424e+00,
        1.32653e+00,  1.37008e+00,  1.46768e+00,  1.37669e+00,
        1.48045e+00,  1.42532e+00,  1.53353e+00,  1.65843e+00,
        1.61700e+00,  1.50176e+00,  1.61346e+00,  1.59878e+00,
        1.72295e+00,  1.57071e+00,  1.65868e+00,  1.61570e+00,
        1.70525e+00,  1.82307e+00,  1.66712e+00,  1.78903e+00,
        1.65301e+00,  1.86202e+00,  1.70566e+00,  1.86215e+00,
        1.73174e+00,  1.62544e+00,  1.61329e+00,  1.45834e+00,
        1.56535e+00,  1.69653e+00,  1.54077e+00,  1.69841e+00,
        1.49033e+00,  1.62786e+00,  1.50735e+00,  1.66511e+00,
        1.54919e+00,  1.46170e+00,  1.50656e+00,  1.42070e+00,
        1.57528e+00,  1.45355e+00,  1.47089e+00,  1.35960e+00,
        1.47545e+00,  1.51723e+00,  1.39286e+00,  1.28524e+00,
        1.34085e+00,  1.23776e+00,  1.32933e+00,  1.23248e+00,
        1.18952e+00,  1.27822e+00,  1.12289e+00,  1.19641e+00,
        1.11737e+00,  1.00946e+00,  1.05478e+00,  1.00369e+00,
        9.34490e-01,  9.22730e-01,  9.65960e-01,  9.14760e-01,
        8.77650e-01,  8.00630e-01,  7.39740e-01,  7.07350e-01,
        6.83780e-01,  6.88370e-01,  6.11130e-01,  5.49760e-01,
        5.38530e-01,  4.86690e-01,  4.61320e-01,  4.30660e-01,
        4.07840e-01,  3.85930e-01,  3.38830e-01,  3.42830e-01,
        3.03300e-01,  2.73390e-01,  2.50600e-01,  2.38570e-01,
        2.12730e-01,  1.91130e-01,  1.73980e-01,  1.60080e-01,
        1.52690e-01,  1.32410e-01,  1.18750e-01,  1.09360e-01,
        9.84100e-02,  8.85800e-02,  7.57500e-02,  6.80400e-02,
        5.91800e-02,  5.20700e-02,  4.46100e-02,  3.97600e-02,
        3.45100e-02,  3.04600e-02,  2.75800e-02,  2.48800e-02,
        2.23000e-02,  1.98700e-02,  1.77200e-02,  1.56800e-02,
        1.40200e-02,  1.24100e-02,  1.09100e-02,  7.02000e-03,
        8.37000e-03,  7.59000e-03,  7.04000e-03,  6.77000e-03,
        6.67000e-03,  6.76000e-03,  7.06000e-03,  7.58000e-03,
        8.39000e-03,  7.05000e-03,  1.08800e-02,  1.23100e-02,
        1.37600e-02,  1.50200e-02,  1.60000e-02,  1.55500e-02,
        1.20100e-02,  4.39000e-03, -8.68000e-03, -1.64700e-02,
       -1.89300e-02, -1.82500e-02, -1.77400e-02, -1.37200e-02,
       -1.00400e-02, -4.71000e-03, -1.20000e-04,  9.61000e-03,
        1.62900e-02,  2.38600e-02,  3.03300e-02,  4.06800e-02,
        5.79400e-02,  6.24800e-02,  7.33600e-02,  8.75100e-02,
        1.06380e-01,  1.29310e-01,  1.38350e-01,  1.58600e-01,
        1.85660e-01,  2.22460e-01,  2.15540e-01,  2.59910e-01,
        2.79570e-01,  2.99660e-01,  3.27540e-01,  3.50290e-01,
        3.99680e-01,  4.08800e-01,  4.36520e-01,  5.10780e-01,
        5.03310e-01,  5.23930e-01,  5.53430e-01,  6.11520e-01,
        6.86020e-01,  7.20500e-01,  7.69120e-01,  7.23310e-01,
        7.32970e-01,  7.99780e-01,  8.48540e-01,  8.01260e-01,
        9.07480e-01,  9.84770e-01,  9.09000e-01,  1.06275e+00,
        9.72030e-01,  1.01346e+00,  1.10854e+00,  1.01540e+00,
        1.11685e+00,  1.05981e+00,  1.16594e+00,  1.28888e+00,
        1.24589e+00,  1.12860e+00,  1.23867e+00,  1.22279e+00,
        1.34556e+00,  1.19172e+00,  1.27852e+00,  1.23481e+00,
        1.32342e+00,  1.44042e+00,  1.28343e+00,  1.40476e+00,
        1.26792e+00,  1.47726e+00,  1.32032e+00,  1.47669e+00,
        1.34592e+00,  1.23938e+00,  1.56246e+00,  1.40649e+00,
        1.51231e+00,  1.64232e+00,  1.48522e+00,  1.64087e+00,
        1.43060e+00,  1.56593e+00,  1.44320e+00,  1.59838e+00,
        1.47949e+00,  1.38881e+00,  1.43064e+00,  1.34172e+00,
        1.49323e+00,  1.36802e+00,  1.38124e+00,  1.27076e+00,
        1.38744e+00,  1.43008e+00,  1.30660e+00,  1.19989e+00,
        1.25645e+00,  1.15433e+00,  1.24689e+00,  1.15106e+00,
        1.10904e+00,  1.19877e+00,  1.04450e+00,  1.11930e+00,
        1.04140e+00,  9.34650e-01,  9.81150e-01,  9.31270e-01,
        8.63300e-01,  8.52790e-01,  8.97290e-01,  8.47390e-01,
        8.11590e-01,  7.35900e-01,  6.76360e-01,  6.45350e-01,
        6.23160e-01,  6.29170e-01,  5.53210e-01,  5.54210e-01,
        5.43060e-01,  4.91290e-01,  4.66000e-01,  4.35410e-01,
        4.12660e-01,  3.90830e-01,  3.43790e-01,  3.47860e-01,
        3.08390e-01,  2.78550e-01,  2.55820e-01,  2.43860e-01,
        2.18070e-01,  1.96530e-01,  1.79440e-01,  1.65590e-01,
        1.58250e-01,  1.38030e-01,  1.24420e-01,  1.15070e-01,
        1.04160e-01,  9.43800e-02,  8.15900e-02,  7.39200e-02,
        6.51000e-02,  5.80300e-02,  5.06000e-02,  4.57900e-02,
        4.05700e-02,  3.65400e-02,  3.36900e-02,  3.10200e-02,
        2.84600e-02,  2.60500e-02,  2.39300e-02,  2.19100e-02,
        2.02600e-02,  1.86600e-02,  1.71800e-02,  1.33000e-02,
        1.46600e-02,  1.38800e-02,  1.33400e-02,  1.30700e-02,
        1.29700e-02};

real radius = 1.4;
real shaftd = 0.1;
real vel_inf = 10;
real lambda_inf = 1.6;
real rot_vel;
//real timestep = 0.01;
real rho = 1.14;
real chord_len = 0.42;
real thickness = 0.15*0.42;
real wings = 3;
real grx = 5;
real gry = 5;
real chord_approx = 0.42;
real md[] = {0.25,0.25,0.25,0.25};
real sd[] = {0.2,0.1};
/* general use functions defined here */

#define sign(val) val==0 ? 0 : val/fabs(val) /* can create a singularity issue */

/* Analytical Solution of Theta Defined here*/

#define theta_inst(w,t) w*t
#define x_new(r,theta) r*cos(theta)
#define y_new(r,theta) r*sin(theta) 
#define x_new_tip(r,theta) r*cos(theta) - (chord_approx/3)*sin(theta)
#define y_new_tip(r,theta) r*sin(theta) + (chord_approx/3)*cos(theta)
#define x_new_mid(r,theta) r*cos(theta) + (chord_approx/3)*sin(theta)
#define y_new_mid(r,theta) r*sin(theta) - (chord_approx/3)*cos(theta)
#define x_new_end(r,theta) r*cos(theta) + (2*chord_approx/3)*sin(theta)
#define y_new_end(r,theta) r*sin(theta) - (2*chord_approx/3)*cos(theta)

#define beta(theta, u, v) theta + atan(v/u)      /* can create a singularity issue */
#define alpha(lambda, beta) (180/3.14)*acos((lambda - sin(beta))/sqrt(pow((cos(beta)),2) + pow((lambda - sin(beta)),2)))   
#define v_net(u,v) sqrt(u*u + v*v)
#define v_rel(v_net, lambda, beta) v_net*pow((pow((lambda-sin(beta)),2)+pow((cos(beta)),2)),0.5)
#define lambda(v_net) radius*rot_vel/v_net
#define w(lbd) vel_inf*lbd/radius

#define gamma_l  1.4 - 6*(0.06 - (thickness/chord_len))
#define gamma_d  1 - 2.5*(0.06 - (thickness/chord_len)) 
#define alpha_dot(alpha_n, alpha_l, timestep) (0.01744)*((alpha_n - alpha_l)/timestep)  /* can create a singularity issue */
#define kappa(alpha_dot) 0.75 + 0.25*sign(alpha_dot)
#define alpha_ml(alpha,kappa,alpha_dot, relative_vel) (180/3.14)*(alpha*(3.14/180) - gamma_l*kappa*sqrt(fabs(0.5*chord_len*alpha_dot/relative_vel))*sign(alpha_dot))  /* can create a singularity issue */
#define alpha_md(alpha,kappa,alpha_dot, relative_vel) (180/3.14)*(alpha*(3.14/180) - gamma_d*kappa*sqrt(fabs(0.5*chord_len*alpha_dot/relative_vel))*sign(alpha_dot))  /* can create a singularity issue */

#define Cl(alpha) cl_val[alpha] //replace alpha by alpha_ml
#define Cd(alpha) cd_val[alpha] //replacce alpha by alpha_md

#define Ct(C_lm,C_dm,alpha)  fabs(C_lm)*sin(alpha*(3.14/180))-fabs(C_dm)*cos(alpha*(3.14/180))
#define Cn(C_lm,C_dm,alpha)  fabs(C_lm)*cos(alpha*(3.14/180))+fabs(C_dm)*sin(alpha*(3.14/180))

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
#define compxrel(x_cent,x_act) (fabs(x_act - x_cent) > 0.5*(pow(grx,-1)) && fabs(x_act - x_cent) < 1.5*(pow(grx,-1)))
#define compyrel(y_cent,y_act) (fabs(y_act - y_cent) > 0.5*(pow(gry,-1)) && fabs(y_act - y_cent) < 1.5*(pow(gry,-1)))

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

/*
    if (x_vel > 10)
    {
        u = 10;
    }
    else
    {
        u = x_vel;
    }

    if (y_vel > 10)
    {
        v = 10;
    }
    else
    {
        v = x_vel;
    }
*/    
    real u = x_vel;
    real v = y_vel;

    real vnet_act = v_net(u,v);


    real x_vel_last = C_U_M1(c,t);
    real y_vel_last = C_V_M1(c,t);
        
/*    if (time == 0)
    {
        x_vel_last = 0;
        y_vel_last = 0;
    }  */   
     
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
        source = 0.5*rho*u*u*shaftd*grx*gry/4;
	// printf("Center force value added : %g\n", source );
    }

    for (j = 0; j < wings ; j++)
    {
        if (compx(x[0], act_x[j]) && compy(x[1], act_y[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);
         
            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);
 
            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
         
            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];
         
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);

            if(j == 0)
            {
                FILE *data;
                data=fopen("data.txt","a");
                fprintf(data,"%d,%g,%d,%d,%d,%d,%g,%g,%g\n", j,(180/3.14)*act_theta[j],alpha_act,alpha_act_last,alpha_ml_act,alpha_md_act,v_rel_act, Ct_act, Ft_act);
                fclose(data);
            }
                if (act_theta[j] > 6.28*floor(act_theta[j]/6.28) - 3.14/2  && act_theta[j] < 6.28*floor(act_theta[j]/6.28) + 3.14/2) 
                    {
                        source = sd[0]*md[1]*(-Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j])) + source;
                    }
                else
                    {
                        source = sd[0]*md[1]*(Fn_act * cos (act_theta[j]) + Ft_act * sin (act_theta[j]))+ source; 
                    }
            

                

//        // printf("%g\n", source);    
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
        // printf("force-x : %g\n", -source);
        // printf("\n");    
        dS[eqn] = 0;
        }
        if (compxrel(x[0], act_x[j]) && compyrel(x[1], act_y[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);
         
            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);
 
            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
         
            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];
         
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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);
         
            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);
 
            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
         
            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];
         
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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);
         
            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);
 
            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
         
            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];
         
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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);
         
            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);
 
            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
         
            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];
         
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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);
         
            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);
 
            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
         
            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];
         
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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);
         
            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);
 
            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
         
            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];
         
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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);
         
            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);
 
            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
         
            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];
         
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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);
         
            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);
 
            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
         
            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];
         
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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);
         
            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);
 
            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
         
            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];
         
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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);
         
            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);
 
            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
         
            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];
         
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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);
         
            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);
 
            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
         
            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];
         
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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);
         
            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);
 
            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
         
            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];
         
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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);
         
            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);
 
            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
         
            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];
         
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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);
         
            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);
 
            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
         
            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];
         
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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);
         
            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);
 
            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
         
            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];
         
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

    if(j == 0 && source != 0)
    {
        FILE *source_data_x;
        source_data_x=fopen("source_data_x.txt","a");
        fprintf(source_data_x,"%g,%g\n",(180/3.14)*act_theta[j],source);
        fclose(source_data_x);
    }

    }


/*
    if (fabs(source) > 1000000)
    {
        source = sign(source)*1000000;
        // printf("Source val limited to 1e+6\n");
    } 
*/
//    // printf("%g %g\n", compx(x[0], act_x[j]), compy(x[1], act_y[j]));


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

    real u = x_vel;
    real v = y_vel;
/*
    if (x_vel > 10)
    {
        u = 10;
    }
    else
    {
        u = x_vel;
    }

    if (y_vel > 10)
    {
        v = 10;
    }
    else
    {
        v = x_vel;
    }
*/

    real vnet_act = v_net(u,v);
    real x_vel_last = C_U_M1(c,t);
    real y_vel_last = C_V_M1(c,t);
        
/*    if (time == 0)
    {
        x_vel_last = 0;
        y_vel_last = 0;
    }  */   
     
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
        act_x_end[i] = x_new_end(radius, act_theta[i]);
        act_y_end[i] = y_new_end(radius, act_theta[i]);
    }

    C_CENTROID(x,c,t);

    for (j = 0; j < wings ; j++)
    {
        if (compx(x[0], act_x[j]) && compy(x[1], act_y[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);

            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);

            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);

            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];

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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);

            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);

            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);

            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];

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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);

            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);

            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);

            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];

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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);

            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);

            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);

            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];

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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);

            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);

            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);

            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];

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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);

            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);

            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);

            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];

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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);

            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);

            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);

            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];

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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);

            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);

            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);

            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];

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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);

            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);

            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);

            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];

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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);

            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);

            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);

            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];

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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);

            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);

            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);

            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];

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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);

            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);

            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);

            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];

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

        if (compx(x[0], act_x_end[j]) && compy(x[1], act_y_end[j]))
        {
            real lbd_act = lambda(vnet_act);
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);

            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);

            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);

            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];

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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);

            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);

            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);

            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];

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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);

            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);

            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);

            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];

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
            real lbd_act_last = lambda(v_net_last);

            real beta_act = beta(act_theta[j], u, v);
            real beta_act_last = beta(last_theta[j], u_last, v_last);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            int alpha_act = alpha(lbd_act, beta_act);                
            int alpha_act_last = alpha(lbd_act_last, beta_act_last);

            real alpha_dot_act = alpha_dot(alpha_act, alpha_act_last, timestep);

            real kappa_act = kappa(alpha_dot_act);
            int alpha_ml_act = alpha_ml(alpha_act, kappa_act, alpha_dot_act, v_rel_act);         
            int alpha_md_act = alpha_md(alpha_act, kappa_act, alpha_dot_act, v_rel_act);

            real Cl_act = Cl(alpha_ml_act);
            real Cd_act = Cd(alpha_md_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);

            //real Cx_act = Cx[alpha_act];
            //real Cy_act = Cy[alpha_act];

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

    if(j == 0 && source != 0)
    {
        FILE *source_data_y;
        source_data_y=fopen("source_data_y.txt","a");
        fprintf(source_data_y,"%g,%g\n",(180/3.14)*act_theta[j],source);
        fclose(source_data_y);
    }

    }

/*
    if (fabs(source) >1000000)
    {
        source = sign(source)*1000000;
        // printf("Source val limited to 1e+6\n");
    } 
*/
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