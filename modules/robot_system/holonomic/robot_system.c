#include <holonomic/robot_system.h>
#include <math.h>
#include <aversive/error.h>
#include <stdlib.h>

#define DEG(x) ((x) * (180.0 / M_PI))
#define RAD(x) ((x) * (M_PI / 180.0))

void rsh_init(struct robot_system_holonomic *rs) {

    rs->pos = NULL;
    rs->motors[0] = NULL;
    rs->motors[1] = NULL;
    rs->motors[2] = NULL;
    rs->speed = 0.;
    rs->rotation_speed = 0.;
    rs->direction = 0.;
}

void rsh_set_cs(struct robot_system_holonomic *rs, int index, struct cs *cs) {
    rs->motors[index] = cs;
} 


void rsh_set_direction(struct robot_system_holonomic *rs, double theta_v) {
    /** @todo should we do a modulo pi here ? */
    rs->direction = theta_v;
}

void rsh_set_direction_int(void * data, int32_t theta_v_deg){
    struct robot_system_holonomic * rs = data;
    rs->direction = RAD(theta_v_deg);
}

void rsh_set_speed(void *data, int32_t speed) {
    struct robot_system_holonomic *rs = data;
    rs->speed = speed;
}

void rsh_set_rotation_speed(void *data, int32_t omega_d) {
    struct robot_system_holonomic *rs = data;
    rs->rotation_speed = RAD(omega_d);
}

void rsh_update(struct robot_system_holonomic *rs) {
    if(rs->pos == NULL)
        return;

    /** @todo add real physical parameters. */
    double theta_r = holonomic_position_get_a_rad_double(rs->pos);
    int i;

    for(i = 0; i < 3; i++) {
        /** omega = omega_r + omega_t = speed of the wheel i */
        
        /** < the part for translation */
        const double omega_t = - rs->speed / rs->pos->geometry.wheel_radius[i] * 
                     cos(theta_r - rs->direction + rs->pos->geometry.beta[i] - M_PI_2);
        /** < the part for rotation */
        const double omega_r = - rs->rotation_speed * rs->pos->geometry.wheel_distance[i] / rs->pos->geometry.wheel_radius[i];
        
        /** @bug @todo le facteur 10 : on ne le met pas et on a des vitesse 10 x trop haute */
        cs_set_consign(rs->motors[i],((omega_t + omega_r) / M_2_PI * (double)rs->pos->geometry.encoder_resolution / (double)rs->pos->update_frequency)/10);
        //DEBUG(E_ROBOT_SYSTEM, "wheel %d : omega_t=%.1f omega_r=%.1f", i, omega_t, omega_r);
    }
}


void rsh_set_position_manager(struct robot_system_holonomic *rs,
                              struct holonomic_robot_position *pos){
    rs->pos = pos;
}
