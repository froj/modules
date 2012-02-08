/*
 *  Copyright Droids Corporation (2007)
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *  Revision : $Id: blocking_detection_manager.h,v 1.1.2.6 2008-05-09 08:25:10 zer0 Exp $
 *
 *  Olivier MATZ <zer0@droids-corp.org>
 */

/* blocking detection manager */

#ifndef BLOCKING_DETECTION_MANAGER_H_
#define BLOCKING_DETECTION_MANAGER_H_

#include <aversive.h>
#include <control_system_manager.h>


/* detect blocking based on motor current.
 * triggers the blocking if:
 *   - the current in the motor is a above a threshold
 *     during n tests
 *   - the speed is below the threshold (if specified)
 *
 * We suppose that i = k1.V - k2.w
 * (V is the voltage applied on the motor, and w the current speed
 * of the motor)
 */

struct blocking_detection {
	struct cs *cs;
	uint16_t cpt_thres;
	uint16_t cpt;
	uint16_t err_thres;
	int32_t old_feedback;
	int32_t old_consign;

};

/** init module, give the cs as parameter */
void bd_init(struct blocking_detection *bd, struct cs *cs);


void bd_set_thresholds(struct blocking_detection *bd, uint16_t err_thres, uint16_t cpt_thres);

/** reset the blocking */
void bd_reset(struct blocking_detection *bd);

/** function to be called periodically */
void bd_manage(struct blocking_detection *bd);

/** get value of blocking detection */
uint8_t bd_get(struct blocking_detection *bd);

#endif
