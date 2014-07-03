
/*
 * libinter - Interpolation methods library
 *
 * Copyright (c) 2013-2014 FOXEL SA - http://foxel.ch
 * Please read <http://foxel.ch/license> for more information.
 *
 *
 * Author(s):
 *
 * Nils Hamel <nils.hamel@foxel.ch>
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Additional Terms:
 *
 * You are required to preserve legal notices and author attributions in
 * that material or in the Appropriate Legal Notices displayed by works
 * containing it.
 *
 * You are required to attribute the work as explained in the "Usage and
 * Attribution" section of <http://foxel.ch/license>.
 */


/* 
    Source - Includes
 */

    # include "inter-bilinear.h"

/*
    Source - Fast bilinear interpolation method
 */

    inter_C8_t inter_bilinearf(

        inter_C8_t *  bm, 
        inter_Index_t bm_w,
        inter_Index_t bm_h,
        inter_Index_t bm_d, 
        inter_Index_t bm_c,
        inter_Real_t  bm_x,
        inter_Real_t  bm_y

    ) {

        /* Interpolation vector */
        static inter_Real_t cv[4];
        static inter_Real_t cc[4];

        /* Interpolation values */
        static inter_Real_t rx = 0.0;
        static inter_Real_t ry = 0.0;

        /* Interpolation variables */
        static inter_Index_t px = 0;
        static inter_Index_t py = 0;

        /* Interpolated value */
        static inter_Real_t iv = 0.0;

        /* Compute relative grid parameters */
        px = trunc( bm_x ); rx = bm_x - px;
        py = trunc( bm_y ); ry = bm_y - py;

        /* Compute interpolation vector */
        cv[0] = * ( bm + bm_d * ( bm_w * ( py    ) + ( px     ) ) + bm_c );
        cv[1] = * ( bm + bm_d * ( bm_w * ( py ++ ) + ( px + 1 ) ) + bm_c );
        cv[2] = * ( bm + bm_d * ( bm_w * ( py    ) + ( px     ) ) + bm_c );
        cv[3] = * ( bm + bm_d * ( bm_w * ( py    ) + ( px + 1 ) ) + bm_c );

        /* Compute interpolation matrix product */
        cc[0] = + cv[ 0];
        cc[1] = - cv[ 0] + cv[ 2];
        cc[2] = - cv[ 0] + cv[ 1];
        cc[3] = + cv[ 0] - cv[ 1] - cv[ 2] + cv[ 3];

        /* Compute interpolated value */
        iv = cc[0] + cc[1] * ry + cc[2] * rx + cc[3] * rx * ry;

        /* Verify interpolated value */
        iv = ( iv <   0.0 ) ?   0.0 : iv; 
        iv = ( iv > 255.0 ) ? 255.0 : iv;

        /* Return interpolated value */
        return( iv );

    }

