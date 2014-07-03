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

    # include "inter-bicubic.h"

/*
    Source - Fast bicubic interpolation method
 */

    inter_C8_t inter_bicubicf(

        inter_C8_t *  bm, 
        inter_Index_t bm_w,
        inter_Index_t bm_h,
        inter_Index_t bm_d, 
        inter_Index_t bm_c,
        inter_Real_t  bm_x,
        inter_Real_t  bm_y

    ) {

        /* Interpolation vector */
        static inter_Real_t cv[16];
        static inter_Real_t cc[16];

        /* Interpolation values */
        static inter_Real_t rx1 = 0.0;
        static inter_Real_t ry1 = 0.0;
        static inter_Real_t rx2 = 0.0;
        static inter_Real_t ry2 = 0.0;
        static inter_Real_t rx3 = 0.0;
        static inter_Real_t ry3 = 0.0;

        /* Interpolation variables */
        static inter_Index_t px = 0;
        static inter_Index_t py = 0;

        /* Interpolated value */
        static inter_Index_t pxm1 = 0;
        static inter_Index_t pxp1 = 0;
        static inter_Index_t pxp2 = 0;
        static inter_Index_t pym1 = 0;
        static inter_Index_t pyp1 = 0;
        static inter_Index_t pyp2 = 0;

        /* Interpolated value */
        static inter_Real_t iv = 0.0;

        /* Compute relative grid parameters */
        px = trunc( bm_x );
        py = trunc( bm_y );

        /* Compute sampling nodes */
        pxm1 = px - 1; pxm1 = ( ( pxm1 <  0    ) ? 0 : pxm1 );
        pxp1 = px + 1;
        pxp2 = px + 2; pxp2 = ( ( pxp2 >= bm_w ) ? bm_w - 1 : pxp2 );
        pym1 = py - 1; pym1 = ( ( pym1 <  0    ) ? 0 : pym1 );
        pyp1 = py + 1;
        pyp2 = py + 2; pyp2 = ( ( pyp2 >= bm_h ) ? bm_h - 1 : pyp2 );

        /* Compute interpolation vector */
        cv[ 0] = * ( bm + bm_d * ( bm_w * pym1 + pxm1 ) + bm_c );
        cv[ 1] = * ( bm + bm_d * ( bm_w * pym1 + px   ) + bm_c );
        cv[ 2] = * ( bm + bm_d * ( bm_w * pym1 + pxp1 ) + bm_c );
        cv[ 3] = * ( bm + bm_d * ( bm_w * pym1 + pxp2 ) + bm_c );
        cv[ 4] = * ( bm + bm_d * ( bm_w * py   + pxm1 ) + bm_c );
        cv[ 5] = * ( bm + bm_d * ( bm_w * py   + px   ) + bm_c );
        cv[ 6] = * ( bm + bm_d * ( bm_w * py   + pxp1 ) + bm_c );
        cv[ 7] = * ( bm + bm_d * ( bm_w * py   + pxp2 ) + bm_c );
        cv[ 8] = * ( bm + bm_d * ( bm_w * pyp1 + pxm1 ) + bm_c );
        cv[ 9] = * ( bm + bm_d * ( bm_w * pyp1 + px   ) + bm_c );
        cv[10] = * ( bm + bm_d * ( bm_w * pyp1 + pxp1 ) + bm_c );
        cv[11] = * ( bm + bm_d * ( bm_w * pyp1 + pxp2 ) + bm_c );
        cv[12] = * ( bm + bm_d * ( bm_w * pyp2 + pxm1 ) + bm_c );
        cv[13] = * ( bm + bm_d * ( bm_w * pyp2 + px   ) + bm_c );
        cv[14] = * ( bm + bm_d * ( bm_w * pyp2 + pxp1 ) + bm_c );
        cv[15] = * ( bm + bm_d * ( bm_w * pyp2 + pxp2 ) + bm_c );

        /* Compute interpolation matrix product */
        cc[ 0] = (   36.0 /   36.0 ) * cv[ 0]; 
        cc[ 1] = (  -66.0 /   36.0 ) * cv[ 0] + (  108.0 /   36.0 ) * cv[ 4] + 
                 (  -54.0 /   36.0 ) * cv[ 8] + (   12.0 /   36.0 ) * cv[12]; 
        cc[ 2] = (   36.0 /   36.0 ) * cv[ 0] + (  -90.0 /   36.0 ) * cv[ 4] + 
                 (   72.0 /   36.0 ) * cv[ 8] + (  -18.0 /   36.0 ) * cv[12]; 
        cc[ 3] = (   -6.0 /   36.0 ) * cv[ 0] + (   18.0 /   36.0 ) * cv[ 4] + 
                 (  -18.0 /   36.0 ) * cv[ 8] + (    6.0 /   36.0 ) * cv[12];
        cc[ 4] = (  -66.0 /   36.0 ) * cv[ 0] + (  108.0 /   36.0 ) * cv[ 1] + 
                 (  -54.0 /   36.0 ) * cv[ 2] + (   12.0 /   36.0 ) * cv[ 3] + 
                 (   -0.0 /   36.0 ) * cv[ 4] + (    0.0 /   36.0 ) * cv[ 7] + 
                 (    0.0 /   36.0 ) * cv[12] + (   -0.0 /   36.0 ) * cv[15];
        cc[ 5] = (  121.0 /   36.0 ) * cv[ 0] + ( -198.0 /   36.0 ) * cv[ 1] + 
                 (   99.0 /   36.0 ) * cv[ 2] + (  -22.0 /   36.0 ) * cv[ 3] + 
                 ( -198.0 /   36.0 ) * cv[ 4] + (  324.0 /   36.0 ) * cv[ 5] + 
                 ( -162.0 /   36.0 ) * cv[ 6] + (   36.0 /   36.0 ) * cv[ 7] + 
                 (   99.0 /   36.0 ) * cv[ 8] + ( -162.0 /   36.0 ) * cv[ 9] + 
                 (   81.0 /   36.0 ) * cv[10] + (  -18.0 /   36.0 ) * cv[11] + 
                 (  -22.0 /   36.0 ) * cv[12] + (   36.0 /   36.0 ) * cv[13] + 
                 (  -18.0 /   36.0 ) * cv[14] + (    4.0 /   36.0 ) * cv[15];
        cc[ 6] = (  -66.0 /   36.0 ) * cv[ 0] + (  108.0 /   36.0 ) * cv[ 1] + 
                 (  -54.0 /   36.0 ) * cv[ 2] + (   12.0 /   36.0 ) * cv[ 3] + 
                 (  165.0 /   36.0 ) * cv[ 4] + ( -270.0 /   36.0 ) * cv[ 5] + 
                 (  135.0 /   36.0 ) * cv[ 6] + (  -30.0 /   36.0 ) * cv[ 7] + 
                 ( -132.0 /   36.0 ) * cv[ 8] + (  216.0 /   36.0 ) * cv[ 9] + 
                 ( -108.0 /   36.0 ) * cv[10] + (   24.0 /   36.0 ) * cv[11] + 
                 (   33.0 /   36.0 ) * cv[12] + (  -54.0 /   36.0 ) * cv[13] + 
                 (   27.0 /   36.0 ) * cv[14] + (   -6.0 /   36.0 ) * cv[15];
        cc[ 7] = (   11.0 /   36.0 ) * cv[ 0] + (  -18.0 /   36.0 ) * cv[ 1] + 
                 (    9.0 /   36.0 ) * cv[ 2] + (   -2.0 /   36.0 ) * cv[ 3] + 
                 (  -33.0 /   36.0 ) * cv[ 4] + (   54.0 /   36.0 ) * cv[ 5] + 
                 (  -27.0 /   36.0 ) * cv[ 6] + (    6.0 /   36.0 ) * cv[ 7] + 
                 (   33.0 /   36.0 ) * cv[ 8] + (  -54.0 /   36.0 ) * cv[ 9] + 
                 (   27.0 /   36.0 ) * cv[10] + (   -6.0 /   36.0 ) * cv[11] + 
                 (  -11.0 /   36.0 ) * cv[12] + (   18.0 /   36.0 ) * cv[13] + 
                 (   -9.0 /   36.0 ) * cv[14] + (    2.0 /   36.0 ) * cv[15];
        cc[ 8] = (   36.0 /   36.0 ) * cv[ 0] + (  -90.0 /   36.0 ) * cv[ 1] + 
                 (   72.0 /   36.0 ) * cv[ 2] + (  -18.0 /   36.0 ) * cv[ 3];
        cc[ 9] = (  -66.0 /   36.0 ) * cv[ 0] + (  165.0 /   36.0 ) * cv[ 1] + 
                 ( -132.0 /   36.0 ) * cv[ 2] + (   33.0 /   36.0 ) * cv[ 3] + 
                 (  108.0 /   36.0 ) * cv[ 4] + ( -270.0 /   36.0 ) * cv[ 5] + 
                 (  216.0 /   36.0 ) * cv[ 6] + (  -54.0 /   36.0 ) * cv[ 7] + 
                 (  -54.0 /   36.0 ) * cv[ 8] + (  135.0 /   36.0 ) * cv[ 9] + 
                 ( -108.0 /   36.0 ) * cv[10] + (   27.0 /   36.0 ) * cv[11] + 
                 (   12.0 /   36.0 ) * cv[12] + (  -30.0 /   36.0 ) * cv[13] + 
                 (   24.0 /   36.0 ) * cv[14] + (   -6.0 /   36.0 ) * cv[15];
        cc[10] = (   36.0 /   36.0 ) * cv[ 0] + (  -90.0 /   36.0 ) * cv[ 1] + 
                 (   72.0 /   36.0 ) * cv[ 2] + (  -18.0 /   36.0 ) * cv[ 3] + 
                 (  -90.0 /   36.0 ) * cv[ 4] + (  225.0 /   36.0 ) * cv[ 5] + 
                 ( -180.0 /   36.0 ) * cv[ 6] + (   45.0 /   36.0 ) * cv[ 7] + 
                 (   72.0 /   36.0 ) * cv[ 8] + ( -180.0 /   36.0 ) * cv[ 9] + 
                 (  144.0 /   36.0 ) * cv[10] + (  -36.0 /   36.0 ) * cv[11] + 
                 (  -18.0 /   36.0 ) * cv[12] + (   45.0 /   36.0 ) * cv[13] + 
                 (  -36.0 /   36.0 ) * cv[14] + (    9.0 /   36.0 ) * cv[15];
        cc[11] = (   -6.0 /   36.0 ) * cv[ 0] + (   15.0 /   36.0 ) * cv[ 1] + 
                 (  -12.0 /   36.0 ) * cv[ 2] + (    3.0 /   36.0 ) * cv[ 3] + 
                 (   18.0 /   36.0 ) * cv[ 4] + (  -45.0 /   36.0 ) * cv[ 5] + 
                 (   36.0 /   36.0 ) * cv[ 6] + (   -9.0 /   36.0 ) * cv[ 7] + 
                 (  -18.0 /   36.0 ) * cv[ 8] + (   45.0 /   36.0 ) * cv[ 9] + 
                 (  -36.0 /   36.0 ) * cv[10] + (    9.0 /   36.0 ) * cv[11] + 
                 (    6.0 /   36.0 ) * cv[12] + (  -15.0 /   36.0 ) * cv[13] + 
                 (   12.0 /   36.0 ) * cv[14] + (   -3.0 /   36.0 ) * cv[15];
        cc[12] = (   -6.0 /   36.0 ) * cv[ 0] + (   18.0 /   36.0 ) * cv[ 1] + 
                 (  -18.0 /   36.0 ) * cv[ 2] + (    6.0 /   36.0 ) * cv[ 3];
        cc[13] = (   11.0 /   36.0 ) * cv[ 0] + (  -33.0 /   36.0 ) * cv[ 1] + 
                 (   33.0 /   36.0 ) * cv[ 2] + (  -11.0 /   36.0 ) * cv[ 3] + 
                 (  -18.0 /   36.0 ) * cv[ 4] + (   54.0 /   36.0 ) * cv[ 5] + 
                 (  -54.0 /   36.0 ) * cv[ 6] + (   18.0 /   36.0 ) * cv[ 7] + 
                 (    9.0 /   36.0 ) * cv[ 8] + (  -27.0 /   36.0 ) * cv[ 9] + 
                 (   27.0 /   36.0 ) * cv[10] + (   -9.0 /   36.0 ) * cv[11] + 
                 (   -2.0 /   36.0 ) * cv[12] + (    6.0 /   36.0 ) * cv[13] + 
                 (   -6.0 /   36.0 ) * cv[14] + (    2.0 /   36.0 ) * cv[15];
        cc[14] = (   -6.0 /   36.0 ) * cv[ 0] + (   18.0 /   36.0 ) * cv[ 1] + 
                 (  -18.0 /   36.0 ) * cv[ 2] + (    6.0 /   36.0 ) * cv[ 3] + 
                 (   15.0 /   36.0 ) * cv[ 4] + (  -45.0 /   36.0 ) * cv[ 5] + 
                 (   45.0 /   36.0 ) * cv[ 6] + (  -15.0 /   36.0 ) * cv[ 7] + 
                 (  -12.0 /   36.0 ) * cv[ 8] + (   36.0 /   36.0 ) * cv[ 9] + 
                 (  -36.0 /   36.0 ) * cv[10] + (   12.0 /   36.0 ) * cv[11] + 
                 (    3.0 /   36.0 ) * cv[12] + (   -9.0 /   36.0 ) * cv[13] + 
                 (    9.0 /   36.0 ) * cv[14] + (   -3.0 /   36.0 ) * cv[15];
        cc[15] = (    1.0 /   36.0 ) * cv[ 0] + (   -3.0 /   36.0 ) * cv[ 1] + 
                 (    3.0 /   36.0 ) * cv[ 2] + (   -1.0 /   36.0 ) * cv[ 3] + 
                 (   -3.0 /   36.0 ) * cv[ 4] + (    9.0 /   36.0 ) * cv[ 5] + 
                 (   -9.0 /   36.0 ) * cv[ 6] + (    3.0 /   36.0 ) * cv[ 7] + 
                 (    3.0 /   36.0 ) * cv[ 8] + (   -9.0 /   36.0 ) * cv[ 9] + 
                 (    9.0 /   36.0 ) * cv[10] + (   -3.0 /   36.0 ) * cv[11] + 
                 (   -1.0 /   36.0 ) * cv[12] + (    3.0 /   36.0 ) * cv[13] + 
                 (   -3.0 /   36.0 ) * cv[14] + (    1.0 /   36.0 ) * cv[15];

        /* Prepare interpolated value computation */
        rx1 = ( bm_x + 1.0 ) - px; rx2 = rx1 * rx1; rx3 = rx1 * rx2;
        ry1 = ( bm_y + 1.0 ) - py; ry2 = ry1 * ry1; ry3 = ry1 * ry2;

        /* Compute interpolated value */
        iv = cc[ 0]       + cc[ 1] * ry1       + cc[ 2] * ry2       + cc[ 3] * ry3       +
             cc[ 4] * rx1 + cc[ 5] * ry1 * rx1 + cc[ 6] * ry2 * rx1 + cc[ 7] * ry3 * rx1 +
             cc[ 8] * rx2 + cc[ 9] * ry1 * rx2 + cc[10] * ry2 * rx2 + cc[11] * ry3 * rx2 +
             cc[12] * rx3 + cc[13] * ry1 * rx3 + cc[14] * ry2 * rx3 + cc[15] * ry3 * rx3;

        /* Verify interpolated value */
        iv = ( iv <   0.0 ) ?   0.0 : iv; 
        iv = ( iv > 255.0 ) ? 255.0 : iv;

        /* Return interpolated value */
        return( iv );

    }

