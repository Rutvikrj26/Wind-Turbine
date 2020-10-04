#include "udf.h"
 
real cl_val [] = {0.    ,  0.11  ,  0.22  ,  0.33  ,  0.418 ,  0.518 ,  0.6048,
        0.676 ,  0.7189,  0.6969,  0.5122,  0.1642,  0.0749,  0.0967,
        0.1382,  0.1861,  0.2364,  0.2873,  0.3393,  0.3927,  0.4463,
        0.5001,  0.5539,  0.6078,  0.6617,  0.7156,  0.77  ,  0.8277,
        0.8368,  0.8459,  0.855 ,  0.88  ,  0.905 ,  0.93  ,  0.955 ,
        0.98  ,  0.991 ,  1.002 ,  1.013 ,  1.024 ,  1.035 ,  1.038 ,
        1.041 ,  1.044 ,  1.047 ,  1.05  ,  1.044 ,  1.038 ,  1.032 ,
        1.026 ,  1.02  ,  1.007 ,  0.994 ,  0.981 ,  0.968 ,  0.955 ,
        0.939 ,  0.923 ,  0.907 ,  0.891 ,  0.875 ,  0.852 ,  0.829 ,
        0.806 ,  0.783 ,  0.76  ,  0.734 ,  0.708 ,  0.682 ,  0.656 ,
        0.63  ,  0.604 ,  0.578 ,  0.552 ,  0.526 ,  0.5   ,  0.473 ,
        0.446 ,  0.419 ,  0.392 ,  0.365 ,  0.338 ,  0.311 ,  0.284 ,
        0.257 ,  0.23  ,  0.202 ,  0.174 ,  0.146 ,  0.118 ,  0.09  ,
        0.062 ,  0.034 ,  0.006 , -0.022 , -0.05  , -0.077 , -0.104 ,
       -0.131 , -0.158 , -0.185 , -0.212 , -0.239 , -0.266 , -0.293 ,
       -0.32  , -0.346 , -0.372 , -0.398 , -0.424 , -0.45  , -0.475 ,
       -0.5   , -0.525 , -0.55  , -0.575 , -0.594 , -0.613 , -0.632 ,
       -0.651 , -0.67  , -0.688 , -0.706 , -0.724 , -0.742 , -0.76  ,
       -0.778 , -0.796 , -0.814 , -0.832 , -0.85  , -0.866 , -0.882 ,
       -0.898 , -0.914 , -0.93  , -0.94  , -0.95  , -0.96  , -0.97  ,
       -0.98  , -0.964 , -0.948 , -0.932 , -0.916 , -0.9   , -0.874 ,
       -0.848 , -0.822 , -0.796 , -0.77  , -0.75  , -0.73  , -0.71  ,
       -0.69  , -0.67  , -0.663 , -0.656 , -0.649 , -0.642 , -0.635 ,
       -0.644 , -0.653 , -0.662 , -0.671 , -0.68  , -0.714 , -0.748 ,
       -0.782 , -0.816 , -0.85  , -0.812 , -0.774 , -0.736 , -0.698 ,
       -0.66  , -0.528 , -0.396 , -0.264 , -0.132};

real cd_val [] = {0.0147    , 0.0148    , 0.0151    , 0.0156    , 0.0168    ,
       0.0181    , 0.0197    , 0.0214    , 0.0234    , 0.0255    ,
       0.0277    , 0.076     , 0.123     , 0.14      , 0.158     ,
       0.177     , 0.196     , 0.217     , 0.238     , 0.26      ,
       0.282     , 0.305     , 0.329     , 0.354     , 0.379     ,
       0.405     , 0.432     , 0.46      , 0.49666667, 0.53333333,
       0.57      , 0.605     , 0.64      , 0.675     , 0.71      ,
       0.745     , 0.78      , 0.815     , 0.85      , 0.885     ,
       0.92      , 0.951     , 0.982     , 1.013     , 1.044     ,
       1.075     , 0.903     , 0.731     , 0.559     , 0.387     ,
       0.215     , 0.441     , 0.667     , 0.893     , 1.119     ,
       1.345     , 1.37      , 1.395     , 1.42      , 1.445     ,
       1.47      , 1.491     , 1.512     , 1.533     , 1.554     ,
       1.575     , 1.593     , 1.611     , 1.629     , 1.647     ,
       1.665     , 1.679     , 1.693     , 1.707     , 1.721     ,
       1.735     , 1.744     , 1.753     , 1.762     , 1.771     ,
       1.78      , 1.784     , 1.788     , 1.792     , 1.796     ,
       1.8       , 1.8       , 1.8       , 1.8       , 1.8       ,
       1.8       , 1.796     , 1.792     , 1.788     , 1.784     ,
       1.78      , 1.774     , 1.768     , 1.762     , 1.756     ,
       1.75      , 1.74      , 1.73      , 1.72      , 1.71      ,
       1.7       , 1.687     , 1.674     , 1.661     , 1.648     ,
       1.635     , 1.619     , 1.603     , 1.587     , 1.571     ,
       1.555     , 1.537     , 1.519     , 1.501     , 1.483     ,
       1.465     , 1.442     , 1.419     , 1.396     , 1.373     ,
       1.35      , 1.325     , 1.3       , 1.275     , 1.25      ,
       1.225     , 1.197     , 1.169     , 1.141     , 1.113     ,
       1.085     , 1.053     , 1.021     , 0.989     , 0.957     ,
       0.925     , 0.891     , 0.857     , 0.823     , 0.789     ,
       0.755     , 0.719     , 0.683     , 0.647     , 0.611     ,
       0.575     , 0.544     , 0.513     , 0.482     , 0.451     ,
       0.42      , 0.4       , 0.38      , 0.36      , 0.34      ,
       0.32      , 0.302     , 0.284     , 0.266     , 0.248     ,
       0.23      , 0.212     , 0.194     , 0.176     , 0.158     ,
       0.14      , 0.123     , 0.106     , 0.089     , 0.072     ,
       0.055     , 0.049     , 0.043     , 0.037     , 0.031     };

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
real grx = 10;
real gry = 10;
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
#define alpha_dot(alpha_n, alpha_l, timestep) ((alpha_n - alpha_l)/timestep)  /* can create a singularity issue */
#define kappa(alpha_dot) 0.75 + 0.25*sign(alpha_dot)
#define alpha_ml(alpha,kappa,alpha_dot, relative_vel) (180/3.14)*(alpha - gamma_l*kappa*sqrt(fabs(0.5*chord_len*alpha_dot/relative_vel))*sign(alpha_dot))  /* can create a singularity issue */
#define alpha_md(alpha,kappa,alpha_dot, relative_vel) (180/3.14)*(alpha - gamma_d*kappa*sqrt(fabs(0.5*chord_len*alpha_dot/relative_vel))*sign(alpha_dot))  /* can create a singularity issue */

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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                

 
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
            real Ct_act = Ct(Cl_act, Cd_act, alpha_act);
            real Cn_act = Cn(Cl_act, Cd_act, alpha_act);
         

         
            real Fn_act = Fn(v_rel_act, Cn_act);                
            real Ft_act = Ft(v_rel_act, Ct_act);

                FILE *data;
                data=fopen("data.txt","a");
                fprintf(data,"%d,%g,%g,%g,%g,%g\n", j,(180/3.14)*act_theta[j],alpha_act,v_rel_act, Ct_act, Ft_act);
                fclose(data);

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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                

 
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                

 
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                
 
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                
 
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                

 
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                

 
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                

 
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                

 
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                

 
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                

 
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                

 
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                

 
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                

 
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                

 
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                

 
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

    real vnet_act = v_net(u,v);


    real x_vel_last = C_U_M1(c,t);
    real y_vel_last = C_V_M1(c,t);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

        if (compx(x[0], act_x_end[j]) && compy(x[1], act_y_end[j]))
        {
            real lbd_act = lambda(vnet_act);

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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

            real beta_act = beta(act_theta[j], u, v);

            real v_rel_act = v_rel(vnet_act, lbd_act, beta_act);                
            
            real alpha_act = alpha(lbd_act, beta_act);                
            real Cl_act = Cl((int) alpha_act);
            real Cd_act = Cd((int) alpha_act);
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