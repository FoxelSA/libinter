/*
 * libinter - Interpolation methods library
 *
 * Copyright (c) 2013-2014 FOXEL SA - http://foxel.ch
 * Please read <http://foxel.ch/license> for more information.
 *
 *
 * Author(s):
 *
 *      Nils Hamel <n.hamel@foxel.ch>
 *
 *
 * This file is part of the FOXEL project <http://foxel.ch>.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Additional Terms:
 *
 *      You are required to preserve legal notices and author attributions in
 *      that material or in the Appropriate Legal Notices displayed by works
 *      containing it.
 *
 *      You are required to attribute the work as explained in the "Usage and
 *      Attribution" section of <http://foxel.ch/license>.
 */


/* 
    Source - Includes
 */

    # include "inter-bilinear.h"

/*
    Source - Fast bilinear interpolation method
 */

    inter_C8_t inter_bilinearf(

        inter_C8_t *  liBytes, 
        inter_Index_t liWidth,
        inter_Index_t liHeight,
        inter_Index_t liLayer, 
        inter_Index_t liChannel,
        inter_Real_t  liX,
        inter_Real_t  liY

    ) {

        /* Interpolation vectors */
        static inter_Real_t liVS[4];
        static inter_Real_t liVC[4];

        /* Optimization variables */
        static inter_Real_t liTX = 0.0;
        static inter_Real_t liTY = 0.0;

        /* Interpolation variables */
        static inter_Index_t liPX = 0;
        static inter_Index_t liPY = 0;

        /* Interpolated variables */
        static inter_Real_t liIV = 0.0;

        /* Compute relatliIVe grid parameters */
        liPX = trunc( liX ); liTX = liX - liPX;
        liPY = trunc( liY ); liTY = liY - liPY;

        /* Compute interpolation vector */
        liVS[0] = * ( liBytes + liLayer * ( liWidth * ( liPY    ) + ( liPX     ) ) + liChannel );
        liVS[1] = * ( liBytes + liLayer * ( liWidth * ( liPY ++ ) + ( liPX + 1 ) ) + liChannel );
        liVS[2] = * ( liBytes + liLayer * ( liWidth * ( liPY    ) + ( liPX     ) ) + liChannel );
        liVS[3] = * ( liBytes + liLayer * ( liWidth * ( liPY    ) + ( liPX + 1 ) ) + liChannel );

        /* Compute interpolation matrix product */
        liVC[0] = + liVS[0];
        liVC[1] = - liVS[0] + liVS[2];
        liVC[2] = - liVS[0] + liVS[1];
        liVC[3] = + liVS[0] - liVS[1] - liVS[2] + liVS[3];

        /* Compute interpolated value */
        liIV = liVC[0] + liVC[1] * liTY + liVC[2] * liTX + liVC[3] * liTX * liTY;

        /* Verify interpolated value */
        liIV = ( liIV <   0.0 ) ?   0.0 : liIV; 
        liIV = ( liIV > 255.0 ) ? 255.0 : liIV;

        /* Return interpolated value */
        return( liIV );

    }

