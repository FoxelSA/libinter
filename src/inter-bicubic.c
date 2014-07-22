/*
 * libinter - Interpolation methods library
 *
 * ColiPYright (c) 2013-2014 FOXEL SA - http://foxel.ch
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
 * You should have receliIVed a coliPY of the GNU Affero General Public License
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

        inter_C8_t *  liBytes, 
        inter_Index_t liWidth,
        inter_Index_t liHeight,
        inter_Index_t liLayer, 
        inter_Index_t liChannel,
        inter_Real_t  liX,
        inter_Real_t  liY

    ) {

        /* Interpolation vectors */
        static inter_Real_t liVS[16];
        static inter_Real_t liVC[16];

        /* Optimization variables */
        static inter_Real_t liTX1 = 0.0;
        static inter_Real_t liTY1 = 0.0;
        static inter_Real_t liTX2 = 0.0;
        static inter_Real_t liTY2 = 0.0;
        static inter_Real_t liTX3 = 0.0;
        static inter_Real_t liTY3 = 0.0;

        /* Interpolation variables */
        static inter_Index_t liPX = 0;
        static inter_Index_t liPY = 0;

        /* Sampling variables */
        static inter_Index_t liPXm1 = 0;
        static inter_Index_t liPXp1 = 0;
        static inter_Index_t liPXp2 = 0;
        static inter_Index_t liPYm1 = 0;
        static inter_Index_t liPYp1 = 0;
        static inter_Index_t liPYp2 = 0;

        /* Interpolated variables */
        static inter_Real_t liIV = 0.0;

        /* Compute relatliIVe grid parameters */
        liPX = trunc( liX );
        liPY = trunc( liY );

        /* Compute sampling nodes */
        liPXm1 = liPX - 1; liPXm1 = ( ( liPXm1 <  0        ) ? 0            : liPXm1 );
        liPXp1 = liPX + 1;
        liPXp2 = liPX + 2; liPXp2 = ( ( liPXp2 >= liWidth  ) ? liWidth - 1  : liPXp2 );
        liPYm1 = liPY - 1; liPYm1 = ( ( liPYm1 <  0        ) ? 0            : liPYm1 );
        liPYp1 = liPY + 1;
        liPYp2 = liPY + 2; liPYp2 = ( ( liPYp2 >= liHeight ) ? liHeight - 1 : liPYp2 );

        /* Compute interpolation vector */
        liVS[ 0] = * ( liBytes + liLayer * ( liWidth * liPYm1 + liPXm1 ) + liChannel );
        liVS[ 1] = * ( liBytes + liLayer * ( liWidth * liPYm1 + liPX   ) + liChannel );
        liVS[ 2] = * ( liBytes + liLayer * ( liWidth * liPYm1 + liPXp1 ) + liChannel );
        liVS[ 3] = * ( liBytes + liLayer * ( liWidth * liPYm1 + liPXp2 ) + liChannel );
        liVS[ 4] = * ( liBytes + liLayer * ( liWidth * liPY   + liPXm1 ) + liChannel );
        liVS[ 5] = * ( liBytes + liLayer * ( liWidth * liPY   + liPX   ) + liChannel );
        liVS[ 6] = * ( liBytes + liLayer * ( liWidth * liPY   + liPXp1 ) + liChannel );
        liVS[ 7] = * ( liBytes + liLayer * ( liWidth * liPY   + liPXp2 ) + liChannel );
        liVS[ 8] = * ( liBytes + liLayer * ( liWidth * liPYp1 + liPXm1 ) + liChannel );
        liVS[ 9] = * ( liBytes + liLayer * ( liWidth * liPYp1 + liPX   ) + liChannel );
        liVS[10] = * ( liBytes + liLayer * ( liWidth * liPYp1 + liPXp1 ) + liChannel );
        liVS[11] = * ( liBytes + liLayer * ( liWidth * liPYp1 + liPXp2 ) + liChannel );
        liVS[12] = * ( liBytes + liLayer * ( liWidth * liPYp2 + liPXm1 ) + liChannel );
        liVS[13] = * ( liBytes + liLayer * ( liWidth * liPYp2 + liPX   ) + liChannel );
        liVS[14] = * ( liBytes + liLayer * ( liWidth * liPYp2 + liPXp1 ) + liChannel );
        liVS[15] = * ( liBytes + liLayer * ( liWidth * liPYp2 + liPXp2 ) + liChannel );

        /* Compute interpolation matrix product */
        liVC[ 0] = (   36.0 /   36.0 ) * liVS[ 0]; 
        liVC[ 1] = (  -66.0 /   36.0 ) * liVS[ 0] + (  108.0 /   36.0 ) * liVS[ 4] + 
                   (  -54.0 /   36.0 ) * liVS[ 8] + (   12.0 /   36.0 ) * liVS[12]; 
        liVC[ 2] = (   36.0 /   36.0 ) * liVS[ 0] + (  -90.0 /   36.0 ) * liVS[ 4] + 
                   (   72.0 /   36.0 ) * liVS[ 8] + (  -18.0 /   36.0 ) * liVS[12]; 
        liVC[ 3] = (   -6.0 /   36.0 ) * liVS[ 0] + (   18.0 /   36.0 ) * liVS[ 4] + 
                   (  -18.0 /   36.0 ) * liVS[ 8] + (    6.0 /   36.0 ) * liVS[12];
        liVC[ 4] = (  -66.0 /   36.0 ) * liVS[ 0] + (  108.0 /   36.0 ) * liVS[ 1] + 
                   (  -54.0 /   36.0 ) * liVS[ 2] + (   12.0 /   36.0 ) * liVS[ 3] + 
                   (   -0.0 /   36.0 ) * liVS[ 4] + (    0.0 /   36.0 ) * liVS[ 7] + 
                   (    0.0 /   36.0 ) * liVS[12] + (   -0.0 /   36.0 ) * liVS[15];
        liVC[ 5] = (  121.0 /   36.0 ) * liVS[ 0] + ( -198.0 /   36.0 ) * liVS[ 1] + 
                   (   99.0 /   36.0 ) * liVS[ 2] + (  -22.0 /   36.0 ) * liVS[ 3] + 
                   ( -198.0 /   36.0 ) * liVS[ 4] + (  324.0 /   36.0 ) * liVS[ 5] + 
                   ( -162.0 /   36.0 ) * liVS[ 6] + (   36.0 /   36.0 ) * liVS[ 7] + 
                   (   99.0 /   36.0 ) * liVS[ 8] + ( -162.0 /   36.0 ) * liVS[ 9] + 
                   (   81.0 /   36.0 ) * liVS[10] + (  -18.0 /   36.0 ) * liVS[11] + 
                   (  -22.0 /   36.0 ) * liVS[12] + (   36.0 /   36.0 ) * liVS[13] + 
                   (  -18.0 /   36.0 ) * liVS[14] + (    4.0 /   36.0 ) * liVS[15];
        liVC[ 6] = (  -66.0 /   36.0 ) * liVS[ 0] + (  108.0 /   36.0 ) * liVS[ 1] + 
                   (  -54.0 /   36.0 ) * liVS[ 2] + (   12.0 /   36.0 ) * liVS[ 3] + 
                   (  165.0 /   36.0 ) * liVS[ 4] + ( -270.0 /   36.0 ) * liVS[ 5] + 
                   (  135.0 /   36.0 ) * liVS[ 6] + (  -30.0 /   36.0 ) * liVS[ 7] + 
                   ( -132.0 /   36.0 ) * liVS[ 8] + (  216.0 /   36.0 ) * liVS[ 9] + 
                   ( -108.0 /   36.0 ) * liVS[10] + (   24.0 /   36.0 ) * liVS[11] + 
                   (   33.0 /   36.0 ) * liVS[12] + (  -54.0 /   36.0 ) * liVS[13] + 
                   (   27.0 /   36.0 ) * liVS[14] + (   -6.0 /   36.0 ) * liVS[15];
        liVC[ 7] = (   11.0 /   36.0 ) * liVS[ 0] + (  -18.0 /   36.0 ) * liVS[ 1] + 
                   (    9.0 /   36.0 ) * liVS[ 2] + (   -2.0 /   36.0 ) * liVS[ 3] + 
                   (  -33.0 /   36.0 ) * liVS[ 4] + (   54.0 /   36.0 ) * liVS[ 5] + 
                   (  -27.0 /   36.0 ) * liVS[ 6] + (    6.0 /   36.0 ) * liVS[ 7] + 
                   (   33.0 /   36.0 ) * liVS[ 8] + (  -54.0 /   36.0 ) * liVS[ 9] + 
                   (   27.0 /   36.0 ) * liVS[10] + (   -6.0 /   36.0 ) * liVS[11] + 
                   (  -11.0 /   36.0 ) * liVS[12] + (   18.0 /   36.0 ) * liVS[13] + 
                   (   -9.0 /   36.0 ) * liVS[14] + (    2.0 /   36.0 ) * liVS[15];
        liVC[ 8] = (   36.0 /   36.0 ) * liVS[ 0] + (  -90.0 /   36.0 ) * liVS[ 1] + 
                   (   72.0 /   36.0 ) * liVS[ 2] + (  -18.0 /   36.0 ) * liVS[ 3];
        liVC[ 9] = (  -66.0 /   36.0 ) * liVS[ 0] + (  165.0 /   36.0 ) * liVS[ 1] + 
                   ( -132.0 /   36.0 ) * liVS[ 2] + (   33.0 /   36.0 ) * liVS[ 3] + 
                   (  108.0 /   36.0 ) * liVS[ 4] + ( -270.0 /   36.0 ) * liVS[ 5] + 
                   (  216.0 /   36.0 ) * liVS[ 6] + (  -54.0 /   36.0 ) * liVS[ 7] + 
                   (  -54.0 /   36.0 ) * liVS[ 8] + (  135.0 /   36.0 ) * liVS[ 9] + 
                   ( -108.0 /   36.0 ) * liVS[10] + (   27.0 /   36.0 ) * liVS[11] + 
                   (   12.0 /   36.0 ) * liVS[12] + (  -30.0 /   36.0 ) * liVS[13] + 
                   (   24.0 /   36.0 ) * liVS[14] + (   -6.0 /   36.0 ) * liVS[15];
        liVC[10] = (   36.0 /   36.0 ) * liVS[ 0] + (  -90.0 /   36.0 ) * liVS[ 1] + 
                   (   72.0 /   36.0 ) * liVS[ 2] + (  -18.0 /   36.0 ) * liVS[ 3] + 
                   (  -90.0 /   36.0 ) * liVS[ 4] + (  225.0 /   36.0 ) * liVS[ 5] + 
                   ( -180.0 /   36.0 ) * liVS[ 6] + (   45.0 /   36.0 ) * liVS[ 7] + 
                   (   72.0 /   36.0 ) * liVS[ 8] + ( -180.0 /   36.0 ) * liVS[ 9] + 
                   (  144.0 /   36.0 ) * liVS[10] + (  -36.0 /   36.0 ) * liVS[11] + 
                   (  -18.0 /   36.0 ) * liVS[12] + (   45.0 /   36.0 ) * liVS[13] + 
                   (  -36.0 /   36.0 ) * liVS[14] + (    9.0 /   36.0 ) * liVS[15];
        liVC[11] = (   -6.0 /   36.0 ) * liVS[ 0] + (   15.0 /   36.0 ) * liVS[ 1] + 
                   (  -12.0 /   36.0 ) * liVS[ 2] + (    3.0 /   36.0 ) * liVS[ 3] + 
                   (   18.0 /   36.0 ) * liVS[ 4] + (  -45.0 /   36.0 ) * liVS[ 5] + 
                   (   36.0 /   36.0 ) * liVS[ 6] + (   -9.0 /   36.0 ) * liVS[ 7] + 
                   (  -18.0 /   36.0 ) * liVS[ 8] + (   45.0 /   36.0 ) * liVS[ 9] + 
                   (  -36.0 /   36.0 ) * liVS[10] + (    9.0 /   36.0 ) * liVS[11] + 
                   (    6.0 /   36.0 ) * liVS[12] + (  -15.0 /   36.0 ) * liVS[13] + 
                   (   12.0 /   36.0 ) * liVS[14] + (   -3.0 /   36.0 ) * liVS[15];
        liVC[12] = (   -6.0 /   36.0 ) * liVS[ 0] + (   18.0 /   36.0 ) * liVS[ 1] + 
                   (  -18.0 /   36.0 ) * liVS[ 2] + (    6.0 /   36.0 ) * liVS[ 3];
        liVC[13] = (   11.0 /   36.0 ) * liVS[ 0] + (  -33.0 /   36.0 ) * liVS[ 1] + 
                   (   33.0 /   36.0 ) * liVS[ 2] + (  -11.0 /   36.0 ) * liVS[ 3] + 
                   (  -18.0 /   36.0 ) * liVS[ 4] + (   54.0 /   36.0 ) * liVS[ 5] + 
                   (  -54.0 /   36.0 ) * liVS[ 6] + (   18.0 /   36.0 ) * liVS[ 7] + 
                   (    9.0 /   36.0 ) * liVS[ 8] + (  -27.0 /   36.0 ) * liVS[ 9] + 
                   (   27.0 /   36.0 ) * liVS[10] + (   -9.0 /   36.0 ) * liVS[11] + 
                   (   -2.0 /   36.0 ) * liVS[12] + (    6.0 /   36.0 ) * liVS[13] + 
                   (   -6.0 /   36.0 ) * liVS[14] + (    2.0 /   36.0 ) * liVS[15];
        liVC[14] = (   -6.0 /   36.0 ) * liVS[ 0] + (   18.0 /   36.0 ) * liVS[ 1] + 
                   (  -18.0 /   36.0 ) * liVS[ 2] + (    6.0 /   36.0 ) * liVS[ 3] + 
                   (   15.0 /   36.0 ) * liVS[ 4] + (  -45.0 /   36.0 ) * liVS[ 5] + 
                   (   45.0 /   36.0 ) * liVS[ 6] + (  -15.0 /   36.0 ) * liVS[ 7] + 
                   (  -12.0 /   36.0 ) * liVS[ 8] + (   36.0 /   36.0 ) * liVS[ 9] + 
                   (  -36.0 /   36.0 ) * liVS[10] + (   12.0 /   36.0 ) * liVS[11] + 
                   (    3.0 /   36.0 ) * liVS[12] + (   -9.0 /   36.0 ) * liVS[13] + 
                   (    9.0 /   36.0 ) * liVS[14] + (   -3.0 /   36.0 ) * liVS[15];
        liVC[15] = (    1.0 /   36.0 ) * liVS[ 0] + (   -3.0 /   36.0 ) * liVS[ 1] + 
                   (    3.0 /   36.0 ) * liVS[ 2] + (   -1.0 /   36.0 ) * liVS[ 3] + 
                   (   -3.0 /   36.0 ) * liVS[ 4] + (    9.0 /   36.0 ) * liVS[ 5] + 
                   (   -9.0 /   36.0 ) * liVS[ 6] + (    3.0 /   36.0 ) * liVS[ 7] + 
                   (    3.0 /   36.0 ) * liVS[ 8] + (   -9.0 /   36.0 ) * liVS[ 9] + 
                   (    9.0 /   36.0 ) * liVS[10] + (   -3.0 /   36.0 ) * liVS[11] + 
                   (   -1.0 /   36.0 ) * liVS[12] + (    3.0 /   36.0 ) * liVS[13] + 
                   (   -3.0 /   36.0 ) * liVS[14] + (    1.0 /   36.0 ) * liVS[15];

        /* Prepare interpolated value computation */
        liTX1 = ( liX + 1.0 ) - liPX; liTX2 = liTX1 * liTX1; liTX3 = liTX1 * liTX2;
        liTY1 = ( liY + 1.0 ) - liPY; liTY2 = liTY1 * liTY1; liTY3 = liTY1 * liTY2;

        /* Compute interpolated value */
        liIV = liVC[ 0]         + liVC[ 1] * liTY1         + liVC[ 2] * liTY2         + liVC[ 3] * liTY3         +
               liVC[ 4] * liTX1 + liVC[ 5] * liTY1 * liTX1 + liVC[ 6] * liTY2 * liTX1 + liVC[ 7] * liTY3 * liTX1 +
               liVC[ 8] * liTX2 + liVC[ 9] * liTY1 * liTX2 + liVC[10] * liTY2 * liTX2 + liVC[11] * liTY3 * liTX2 +
               liVC[12] * liTX3 + liVC[13] * liTY1 * liTX3 + liVC[14] * liTY2 * liTX3 + liVC[15] * liTY3 * liTX3;

        /* Verify interpolated value */
        liIV = ( liIV <   0.0 ) ?   0.0 : liIV; 
        liIV = ( liIV > 255.0 ) ? 255.0 : liIV;

        /* Return interpolated value */
        return( liIV );

    }

