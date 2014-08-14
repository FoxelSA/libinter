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

    # include "inter-bipentic.h"

/*
    Source - Fast bipentic interpolation method
 */

    li_C8_t li_bipenticf(

        li_C8_t * liBytes, 
        li_Size_t liWidth,
        li_Size_t liHeight,
        li_Size_t liLayer, 
        li_Size_t liChannel,
        li_Real_t liX,
        li_Real_t liY

    ) {

        /* Interpolation vectors */
        li_Real_t liVS[36] = { li_Real_s( 0.0 ) };
        li_Real_t liVC[36] = { li_Real_s( 0.0 ) };

        /* Optimization variables */
        li_Real_t liTX1 = li_Real_s( 0.0 );
        li_Real_t liTY1 = li_Real_s( 0.0 );
        li_Real_t liTX2 = li_Real_s( 0.0 );
        li_Real_t liTY2 = li_Real_s( 0.0 );
        li_Real_t liTX3 = li_Real_s( 0.0 );
        li_Real_t liTY3 = li_Real_s( 0.0 );
        li_Real_t liTX4 = li_Real_s( 0.0 );
        li_Real_t liTY4 = li_Real_s( 0.0 );
        li_Real_t liTX5 = li_Real_s( 0.0 );
        li_Real_t liTY5 = li_Real_s( 0.0 );

        /* Interpolation variables */
        li_Size_t liPX = li_Size_s( 0 );
        li_Size_t liPY = li_Size_s( 0 );

        /* Sampling variables */
        li_Size_t liPXm2 = li_Size_s( 0 );
        li_Size_t liPXm1 = li_Size_s( 0 );
        li_Size_t liPXp1 = li_Size_s( 0 );
        li_Size_t liPXp2 = li_Size_s( 0 );
        li_Size_t liPXp3 = li_Size_s( 0 );
        li_Size_t liPYm2 = li_Size_s( 0 );
        li_Size_t liPYm1 = li_Size_s( 0 );
        li_Size_t liPYp1 = li_Size_s( 0 );
        li_Size_t liPYp2 = li_Size_s( 0 );
        li_Size_t liPYp3 = li_Size_s( 0 );

        /* Interpolated variables */
        li_Real_t liIV = li_Real_s( 0.0 );

        /* Compute relatliIVe grid parameters */
        liPX = li_Trunc( liX );
        liPY = li_Trunc( liY );

        /* Compute sampling nodes */
        liPXp1 = liPX + li_Size_s( 1 );
        liPYp1 = liPY + li_Size_s( 1 );

        /* Compute sampling nodes */
        liPXm2 = liPX - li_Size_s( 2 ); liPXm2 = ( ( liPXm2 <  li_Size_s( 0 ) ) ? li_Size_s( 0 ) : liPXm2 );
        liPXm1 = liPX - li_Size_s( 1 ); liPXm1 = ( ( liPXm1 <  li_Size_s( 0 ) ) ? li_Size_s( 0 ) : liPXm1 );
        liPYm2 = liPY - li_Size_s( 2 ); liPYm2 = ( ( liPYm2 <  li_Size_s( 0 ) ) ? li_Size_s( 0 ) : liPYm2 );
        liPYm1 = liPY - li_Size_s( 1 ); liPYm1 = ( ( liPYm1 <  li_Size_s( 0 ) ) ? li_Size_s( 0 ) : liPYm1 );

        /* Compute sampling nodes */
        liPXp2 = liPX + li_Size_s( 2 ); liPXp2 = ( ( liPXp2 >= liWidth  ) ? liWidth  - li_Size_s( 1 ) : liPXp2 );
        liPXp3 = liPX + li_Size_s( 3 ); liPXp3 = ( ( liPXp3 >= liWidth  ) ? liWidth  - li_Size_s( 1 ) : liPXp3 );
        liPYp2 = liPY + li_Size_s( 2 ); liPYp2 = ( ( liPYp2 >= liHeight ) ? liHeight - li_Size_s( 1 ) : liPYp2 );
        liPYp3 = liPY + li_Size_s( 3 ); liPYp3 = ( ( liPYp3 >= liHeight ) ? liHeight - li_Size_s( 1 ) : liPYp3 );

        /* Compute interpolation vector */
        liVS[ 0] = * ( liBytes + liLayer * ( liWidth * liPYm2 + liPXm2 ) + liChannel );
        liVS[ 1] = * ( liBytes + liLayer * ( liWidth * liPYm2 + liPXm1 ) + liChannel );
        liVS[ 2] = * ( liBytes + liLayer * ( liWidth * liPYm2 + liPX   ) + liChannel );
        liVS[ 3] = * ( liBytes + liLayer * ( liWidth * liPYm2 + liPXp1 ) + liChannel );
        liVS[ 4] = * ( liBytes + liLayer * ( liWidth * liPYm2 + liPXp2 ) + liChannel );
        liVS[ 5] = * ( liBytes + liLayer * ( liWidth * liPYm2 + liPXp3 ) + liChannel );
        liVS[ 6] = * ( liBytes + liLayer * ( liWidth * liPYm1 + liPXm2 ) + liChannel );
        liVS[ 7] = * ( liBytes + liLayer * ( liWidth * liPYm1 + liPXm1 ) + liChannel );
        liVS[ 8] = * ( liBytes + liLayer * ( liWidth * liPYm1 + liPX   ) + liChannel );
        liVS[ 9] = * ( liBytes + liLayer * ( liWidth * liPYm1 + liPXp1 ) + liChannel );
        liVS[10] = * ( liBytes + liLayer * ( liWidth * liPYm1 + liPXp2 ) + liChannel );
        liVS[11] = * ( liBytes + liLayer * ( liWidth * liPYm1 + liPXp3 ) + liChannel );
        liVS[12] = * ( liBytes + liLayer * ( liWidth * liPY   + liPXm2 ) + liChannel );
        liVS[13] = * ( liBytes + liLayer * ( liWidth * liPY   + liPXm1 ) + liChannel );
        liVS[14] = * ( liBytes + liLayer * ( liWidth * liPY   + liPX   ) + liChannel );
        liVS[15] = * ( liBytes + liLayer * ( liWidth * liPY   + liPXp1 ) + liChannel );
        liVS[16] = * ( liBytes + liLayer * ( liWidth * liPY   + liPXp2 ) + liChannel );
        liVS[17] = * ( liBytes + liLayer * ( liWidth * liPY   + liPXp3 ) + liChannel );
        liVS[18] = * ( liBytes + liLayer * ( liWidth * liPYp1 + liPXm2 ) + liChannel );
        liVS[19] = * ( liBytes + liLayer * ( liWidth * liPYp1 + liPXm1 ) + liChannel );
        liVS[20] = * ( liBytes + liLayer * ( liWidth * liPYp1 + liPX   ) + liChannel );
        liVS[21] = * ( liBytes + liLayer * ( liWidth * liPYp1 + liPXp1 ) + liChannel );
        liVS[22] = * ( liBytes + liLayer * ( liWidth * liPYp1 + liPXp2 ) + liChannel );
        liVS[23] = * ( liBytes + liLayer * ( liWidth * liPYp1 + liPXp3 ) + liChannel );
        liVS[24] = * ( liBytes + liLayer * ( liWidth * liPYp2 + liPXm2 ) + liChannel );
        liVS[25] = * ( liBytes + liLayer * ( liWidth * liPYp2 + liPXm1 ) + liChannel );
        liVS[26] = * ( liBytes + liLayer * ( liWidth * liPYp2 + liPX   ) + liChannel );
        liVS[27] = * ( liBytes + liLayer * ( liWidth * liPYp2 + liPXp1 ) + liChannel );
        liVS[28] = * ( liBytes + liLayer * ( liWidth * liPYp2 + liPXp2 ) + liChannel );
        liVS[29] = * ( liBytes + liLayer * ( liWidth * liPYp2 + liPXp3 ) + liChannel );
        liVS[30] = * ( liBytes + liLayer * ( liWidth * liPYp3 + liPXm2 ) + liChannel );
        liVS[31] = * ( liBytes + liLayer * ( liWidth * liPYp3 + liPXm1 ) + liChannel );
        liVS[32] = * ( liBytes + liLayer * ( liWidth * liPYp3 + liPX   ) + liChannel );
        liVS[33] = * ( liBytes + liLayer * ( liWidth * liPYp3 + liPXp1 ) + liChannel );
        liVS[34] = * ( liBytes + liLayer * ( liWidth * liPYp3 + liPXp2 ) + liChannel );
        liVS[35] = * ( liBytes + liLayer * ( liWidth * liPYp3 + liPXp3 ) + liChannel );

        /* Compute interpolation matrix product */
        liVC[ 0] = ( li_Real_s(   14400.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0];
        liVC[ 1] = ( li_Real_s(  -32880.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(   72000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(  -72000.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(   48000.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(  -18000.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(    2880.0 ) / li_Real_s( 14400.0 ) ) * liVS[30];
        liVC[ 2] = ( li_Real_s(   27000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(  -92400.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(  128400.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(  -93600.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(   36600.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(   -6000.0 ) / li_Real_s( 14400.0 ) ) * liVS[30];
        liVC[ 3] = ( li_Real_s(  -10200.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(   42600.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(  -70800.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(   58800.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(  -24600.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(    4200.0 ) / li_Real_s( 14400.0 ) ) * liVS[30];
        liVC[ 4] = ( li_Real_s(    1800.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(   -8400.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(   15600.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(  -14400.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(    6600.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(   -1200.0 ) / li_Real_s( 14400.0 ) ) * liVS[30];
        liVC[ 5] = ( li_Real_s(    -120.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(     600.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(   -1200.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(    1200.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(    -600.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(     120.0 ) / li_Real_s( 14400.0 ) ) * liVS[30];
        liVC[ 6] = ( li_Real_s(  -32880.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(   72000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(  -72000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(   48000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(  -18000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(    2880.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5];
        liVC[ 7] = ( li_Real_s(   75076.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s( -164400.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(  164400.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s( -109600.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(   41100.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(   -6576.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s( -164400.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(  360000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s( -360000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s(  240000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(  -90000.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(   14400.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(  164400.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s( -360000.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s(  360000.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s( -240000.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s(   90000.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(  -14400.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s( -109600.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(  240000.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s( -240000.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s(  160000.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(  -60000.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(    9600.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(   41100.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(  -90000.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s(   90000.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(  -60000.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(   22500.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(   -3600.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(   -6576.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(   14400.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(  -14400.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(    9600.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(   -3600.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(     576.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[ 8] = ( li_Real_s(  -61650.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(  135000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s( -135000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(   90000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(  -33750.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(    5400.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(  210980.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s( -462000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s(  462000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s( -308000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(  115500.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(  -18480.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s( -293180.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(  642000.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s( -642000.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s(  428000.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s( -160500.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(   25680.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(  213720.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s( -468000.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s(  468000.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s( -312000.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(  117000.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(  -18720.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(  -83570.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(  183000.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s( -183000.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(  122000.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(  -45750.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(    7320.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(   13700.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(  -30000.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(   30000.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(  -20000.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(    7500.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(   -1200.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[ 9] = ( li_Real_s(   23290.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(  -51000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(   51000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(  -34000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(   12750.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(   -2040.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(  -97270.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(  213000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s( -213000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s(  142000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(  -53250.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(    8520.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(  161660.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s( -354000.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s(  354000.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s( -236000.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s(   88500.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(  -14160.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s( -134260.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(  294000.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s( -294000.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s(  196000.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(  -73500.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(   11760.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(   56170.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s( -123000.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s(  123000.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(  -82000.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(   30750.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(   -4920.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(   -9590.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(   21000.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(  -21000.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(   14000.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(   -5250.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(     840.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[10] = ( li_Real_s(   -4110.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(    9000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(   -9000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(    6000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(   -2250.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(     360.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(   19180.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(  -42000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s(   42000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s(  -28000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(   10500.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(   -1680.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(  -35620.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(   78000.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s(  -78000.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s(   52000.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s(  -19500.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(    3120.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(   32880.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(  -72000.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s(   72000.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s(  -48000.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(   18000.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(   -2880.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(  -15070.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(   33000.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s(  -33000.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(   22000.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(   -8250.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(    1320.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(    2740.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(   -6000.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(    6000.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(   -4000.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(    1500.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(    -240.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[11] = ( li_Real_s(     274.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(    -600.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(     600.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(    -400.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(     150.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(     -24.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(   -1370.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(    3000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s(   -3000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s(    2000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(    -750.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(     120.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(    2740.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(   -6000.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s(    6000.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s(   -4000.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s(    1500.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(    -240.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(   -2740.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(    6000.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s(   -6000.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s(    4000.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(   -1500.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(     240.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(    1370.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(   -3000.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s(    3000.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(   -2000.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(     750.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(    -120.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(    -274.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(     600.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(    -600.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(     400.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(    -150.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(      24.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[12] = ( li_Real_s(   27000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(  -92400.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(  128400.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(  -93600.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(   36600.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(   -6000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5];
        liVC[13] = ( li_Real_s(  -61650.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(  210980.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s( -293180.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(  213720.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(  -83570.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(   13700.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(  135000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s( -462000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s(  642000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s( -468000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(  183000.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(  -30000.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s( -135000.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(  462000.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s( -642000.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s(  468000.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s( -183000.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(   30000.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(   90000.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s( -308000.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s(  428000.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s( -312000.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(  122000.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(  -20000.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(  -33750.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(  115500.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s( -160500.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(  117000.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(  -45750.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(    7500.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(    5400.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(  -18480.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(   25680.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(  -18720.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(    7320.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(   -1200.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[14] = ( li_Real_s(   50625.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s( -173250.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(  240750.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s( -175500.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(   68625.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(  -11250.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s( -173250.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(  592900.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s( -823900.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s(  600600.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s( -234850.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(   38500.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(  240750.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s( -823900.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s( 1144900.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s( -834600.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s(  326350.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(  -53500.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s( -175500.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(  600600.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s( -834600.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s(  608400.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s( -237900.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(   39000.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(   68625.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s( -234850.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s(  326350.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s( -237900.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(   93025.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(  -15250.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(  -11250.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(   38500.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(  -53500.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(   39000.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(  -15250.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(    2500.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[15] = ( li_Real_s(  -19125.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(   65450.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(  -90950.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(   66300.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(  -25925.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(    4250.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(   79875.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s( -273350.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s(  379850.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s( -276900.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(  108275.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(  -17750.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s( -132750.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(  454300.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s( -631300.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s(  460200.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s( -179950.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(   29500.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(  110250.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s( -377300.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s(  524300.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s( -382200.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(  149450.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(  -24500.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(  -46125.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(  157850.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s( -219350.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(  159900.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(  -62525.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(   10250.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(    7875.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(  -26950.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(   37450.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(  -27300.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(   10675.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(   -1750.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[16] = ( li_Real_s(    3375.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(  -11550.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(   16050.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(  -11700.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(    4575.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(    -750.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(  -15750.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(   53900.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s(  -74900.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s(   54600.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(  -21350.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(    3500.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(   29250.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s( -100100.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s(  139100.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s( -101400.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s(   39650.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(   -6500.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(  -27000.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(   92400.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s( -128400.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s(   93600.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(  -36600.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(    6000.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(   12375.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(  -42350.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s(   58850.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(  -42900.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(   16775.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(   -2750.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(   -2250.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(    7700.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(  -10700.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(    7800.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(   -3050.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(     500.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[17] = ( li_Real_s(    -225.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(     770.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(   -1070.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(     780.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(    -305.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(      50.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(    1125.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(   -3850.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s(    5350.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s(   -3900.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(    1525.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(    -250.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(   -2250.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(    7700.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s(  -10700.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s(    7800.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s(   -3050.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(     500.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(    2250.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(   -7700.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s(   10700.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s(   -7800.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(    3050.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(    -500.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(   -1125.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(    3850.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s(   -5350.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(    3900.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(   -1525.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(     250.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(     225.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(    -770.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(    1070.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(    -780.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(     305.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(     -50.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[18] = ( li_Real_s(  -10200.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(   42600.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(  -70800.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(   58800.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(  -24600.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(    4200.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5];
        liVC[19] = ( li_Real_s(   23290.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(  -97270.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(  161660.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s( -134260.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(   56170.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(   -9590.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(  -51000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(  213000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s( -354000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s(  294000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s( -123000.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(   21000.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(   51000.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s( -213000.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s(  354000.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s( -294000.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s(  123000.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(  -21000.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(  -34000.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(  142000.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s( -236000.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s(  196000.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(  -82000.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(   14000.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(   12750.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(  -53250.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s(   88500.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(  -73500.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(   30750.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(   -5250.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(   -2040.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(    8520.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(  -14160.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(   11760.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(   -4920.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(     840.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[20] = ( li_Real_s(  -19125.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(   79875.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s( -132750.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(  110250.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(  -46125.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(    7875.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(   65450.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s( -273350.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s(  454300.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s( -377300.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(  157850.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(  -26950.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(  -90950.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(  379850.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s( -631300.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s(  524300.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s( -219350.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(   37450.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(   66300.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s( -276900.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s(  460200.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s( -382200.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(  159900.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(  -27300.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(  -25925.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(  108275.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s( -179950.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(  149450.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(  -62525.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(   10675.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(    4250.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(  -17750.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(   29500.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(  -24500.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(   10250.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(   -1750.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[21] = ( li_Real_s(    7225.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(  -30175.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(   50150.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(  -41650.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(   17425.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(   -2975.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(  -30175.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(  126025.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s( -209450.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s(  173950.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(  -72775.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(   12425.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(   50150.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s( -209450.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s(  348100.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s( -289100.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s(  120950.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(  -20650.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(  -41650.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(  173950.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s( -289100.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s(  240100.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s( -100450.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(   17150.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(   17425.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(  -72775.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s(  120950.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s( -100450.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(   42025.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(   -7175.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(   -2975.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(   12425.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(  -20650.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(   17150.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(   -7175.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(    1225.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[22] = ( li_Real_s(   -1275.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(    5325.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(   -8850.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(    7350.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(   -3075.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(     525.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(    5950.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(  -24850.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s(   41300.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s(  -34300.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(   14350.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(   -2450.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(  -11050.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(   46150.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s(  -76700.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s(   63700.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s(  -26650.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(    4550.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(   10200.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(  -42600.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s(   70800.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s(  -58800.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(   24600.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(   -4200.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(   -4675.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(   19525.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s(  -32450.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(   26950.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(  -11275.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(    1925.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(     850.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(   -3550.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(    5900.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(   -4900.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(    2050.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(    -350.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[23] = ( li_Real_s(      85.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(    -355.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(     590.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(    -490.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(     205.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(     -35.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(    -425.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(    1775.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s(   -2950.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s(    2450.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(   -1025.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(     175.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(     850.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(   -3550.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s(    5900.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s(   -4900.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s(    2050.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(    -350.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(    -850.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(    3550.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s(   -5900.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s(    4900.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(   -2050.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(     350.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(     425.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(   -1775.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s(    2950.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(   -2450.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(    1025.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(    -175.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(     -85.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(     355.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(    -590.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(     490.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(    -205.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(      35.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[24] = ( li_Real_s(    1800.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(   -8400.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(   15600.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(  -14400.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(    6600.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(   -1200.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5];
        liVC[25] = ( li_Real_s(   -4110.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(   19180.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(  -35620.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(   32880.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(  -15070.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(    2740.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(    9000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(  -42000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s(   78000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s(  -72000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(   33000.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(   -6000.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(   -9000.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(   42000.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s(  -78000.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s(   72000.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s(  -33000.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(    6000.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(    6000.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(  -28000.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s(   52000.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s(  -48000.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(   22000.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(   -4000.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(   -2250.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(   10500.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s(  -19500.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(   18000.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(   -8250.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(    1500.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(     360.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(   -1680.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(    3120.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(   -2880.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(    1320.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(    -240.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[26] = ( li_Real_s(    3375.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(  -15750.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(   29250.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(  -27000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(   12375.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(   -2250.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(  -11550.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(   53900.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s( -100100.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s(   92400.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(  -42350.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(    7700.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(   16050.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(  -74900.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s(  139100.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s( -128400.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s(   58850.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(  -10700.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(  -11700.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(   54600.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s( -101400.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s(   93600.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(  -42900.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(    7800.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(    4575.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(  -21350.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s(   39650.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(  -36600.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(   16775.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(   -3050.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(    -750.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(    3500.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(   -6500.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(    6000.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(   -2750.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(     500.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[27] = ( li_Real_s(   -1275.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(    5950.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(  -11050.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(   10200.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(   -4675.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(     850.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(    5325.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(  -24850.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s(   46150.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s(  -42600.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(   19525.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(   -3550.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(   -8850.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(   41300.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s(  -76700.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s(   70800.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s(  -32450.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(    5900.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(    7350.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(  -34300.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s(   63700.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s(  -58800.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(   26950.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(   -4900.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(   -3075.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(   14350.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s(  -26650.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(   24600.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(  -11275.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(    2050.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(     525.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(   -2450.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(    4550.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(   -4200.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(    1925.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(    -350.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[28] = ( li_Real_s(     225.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(   -1050.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(    1950.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(   -1800.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(     825.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(    -150.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(   -1050.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(    4900.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s(   -9100.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s(    8400.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(   -3850.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(     700.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(    1950.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(   -9100.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s(   16900.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s(  -15600.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s(    7150.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(   -1300.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(   -1800.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(    8400.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s(  -15600.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s(   14400.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(   -6600.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(    1200.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(     825.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(   -3850.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s(    7150.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(   -6600.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(    3025.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(    -550.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(    -150.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(     700.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(   -1300.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(    1200.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(    -550.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(     100.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[29] = ( li_Real_s(     -15.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(      70.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(    -130.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(     120.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(     -55.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(      10.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(      75.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(    -350.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s(     650.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s(    -600.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(     275.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(     -50.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(    -150.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(     700.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s(   -1300.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s(    1200.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s(    -550.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(     100.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(     150.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(    -700.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s(    1300.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s(   -1200.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(     550.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(    -100.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(     -75.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(     350.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s(    -650.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(     600.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(    -275.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(      50.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(      15.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(     -70.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(     130.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(    -120.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(      55.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(     -10.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[30] = ( li_Real_s(    -120.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(     600.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(   -1200.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(    1200.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(    -600.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(     120.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5];
        liVC[31] = ( li_Real_s(     274.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(   -1370.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(    2740.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(   -2740.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(    1370.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(    -274.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(    -600.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(    3000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s(   -6000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s(    6000.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(   -3000.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(     600.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(     600.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(   -3000.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s(    6000.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s(   -6000.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s(    3000.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(    -600.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(    -400.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(    2000.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s(   -4000.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s(    4000.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(   -2000.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(     400.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(     150.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(    -750.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s(    1500.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(   -1500.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(     750.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(    -150.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(     -24.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(     120.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(    -240.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(     240.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(    -120.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(      24.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[32] = ( li_Real_s(    -225.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(    1125.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(   -2250.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(    2250.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(   -1125.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(     225.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(     770.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(   -3850.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s(    7700.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s(   -7700.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(    3850.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(    -770.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(   -1070.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(    5350.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s(  -10700.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s(   10700.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s(   -5350.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(    1070.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(     780.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(   -3900.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s(    7800.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s(   -7800.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(    3900.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(    -780.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(    -305.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(    1525.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s(   -3050.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(    3050.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(   -1525.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(     305.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(      50.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(    -250.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(     500.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(    -500.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(     250.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(     -50.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[33] = ( li_Real_s(      85.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(    -425.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(     850.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(    -850.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(     425.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(     -85.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(    -355.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(    1775.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s(   -3550.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s(    3550.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(   -1775.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(     355.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(     590.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(   -2950.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s(    5900.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s(   -5900.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s(    2950.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(    -590.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(    -490.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(    2450.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s(   -4900.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s(    4900.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(   -2450.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(     490.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(     205.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(   -1025.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s(    2050.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(   -2050.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(    1025.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(    -205.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(     -35.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(     175.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(    -350.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(     350.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(    -175.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(      35.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[34] = ( li_Real_s(     -15.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(      75.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(    -150.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(     150.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(     -75.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(      15.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(      70.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(    -350.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s(     700.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s(    -700.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(     350.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(     -70.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(    -130.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(     650.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s(   -1300.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s(    1300.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s(    -650.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(     130.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(     120.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(    -600.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s(    1200.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s(   -1200.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(     600.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(    -120.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(     -55.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(     275.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s(    -550.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(     550.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(    -275.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(      55.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(      10.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(     -50.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(     100.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(    -100.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(      50.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(     -10.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];
        liVC[35] = ( li_Real_s(       1.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 0] +
                   ( li_Real_s(      -5.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 1] +
                   ( li_Real_s(      10.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 2] +
                   ( li_Real_s(     -10.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 3] +
                   ( li_Real_s(       5.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 4] +
                   ( li_Real_s(      -1.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 5] +
                   ( li_Real_s(      -5.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 6] +
                   ( li_Real_s(      25.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 7] +
                   ( li_Real_s(     -50.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 8] +
                   ( li_Real_s(      50.0 ) / li_Real_s( 14400.0 ) ) * liVS[ 9] +
                   ( li_Real_s(     -25.0 ) / li_Real_s( 14400.0 ) ) * liVS[10] +
                   ( li_Real_s(       5.0 ) / li_Real_s( 14400.0 ) ) * liVS[11] +
                   ( li_Real_s(      10.0 ) / li_Real_s( 14400.0 ) ) * liVS[12] +
                   ( li_Real_s(     -50.0 ) / li_Real_s( 14400.0 ) ) * liVS[13] +
                   ( li_Real_s(     100.0 ) / li_Real_s( 14400.0 ) ) * liVS[14] +
                   ( li_Real_s(    -100.0 ) / li_Real_s( 14400.0 ) ) * liVS[15] +
                   ( li_Real_s(      50.0 ) / li_Real_s( 14400.0 ) ) * liVS[16] +
                   ( li_Real_s(     -10.0 ) / li_Real_s( 14400.0 ) ) * liVS[17] +
                   ( li_Real_s(     -10.0 ) / li_Real_s( 14400.0 ) ) * liVS[18] +
                   ( li_Real_s(      50.0 ) / li_Real_s( 14400.0 ) ) * liVS[19] +
                   ( li_Real_s(    -100.0 ) / li_Real_s( 14400.0 ) ) * liVS[20] +
                   ( li_Real_s(     100.0 ) / li_Real_s( 14400.0 ) ) * liVS[21] +
                   ( li_Real_s(     -50.0 ) / li_Real_s( 14400.0 ) ) * liVS[22] +
                   ( li_Real_s(      10.0 ) / li_Real_s( 14400.0 ) ) * liVS[23] +
                   ( li_Real_s(       5.0 ) / li_Real_s( 14400.0 ) ) * liVS[24] +
                   ( li_Real_s(     -25.0 ) / li_Real_s( 14400.0 ) ) * liVS[25] +
                   ( li_Real_s(      50.0 ) / li_Real_s( 14400.0 ) ) * liVS[26] +
                   ( li_Real_s(     -50.0 ) / li_Real_s( 14400.0 ) ) * liVS[27] +
                   ( li_Real_s(      25.0 ) / li_Real_s( 14400.0 ) ) * liVS[28] +
                   ( li_Real_s(      -5.0 ) / li_Real_s( 14400.0 ) ) * liVS[29] +
                   ( li_Real_s(      -1.0 ) / li_Real_s( 14400.0 ) ) * liVS[30] +
                   ( li_Real_s(       5.0 ) / li_Real_s( 14400.0 ) ) * liVS[31] +
                   ( li_Real_s(     -10.0 ) / li_Real_s( 14400.0 ) ) * liVS[32] +
                   ( li_Real_s(      10.0 ) / li_Real_s( 14400.0 ) ) * liVS[33] +
                   ( li_Real_s(      -5.0 ) / li_Real_s( 14400.0 ) ) * liVS[34] +
                   ( li_Real_s(       1.0 ) / li_Real_s( 14400.0 ) ) * liVS[35];

        /* Prepare interpolated value computation */
        liTX1 = ( liX + li_Real_s( 2.0 ) ) - liPX;
        liTX2 = liTX1 * liTX1; 
        liTX3 = liTX1 * liTX2; 
        liTX4 = liTX1 * liTX3; 
        liTX5 = liTX1 * liTX4;
        liTY1 = ( liY + li_Real_s( 2.0 ) ) - liPY; 
        liTY2 = liTY1 * liTY1; 
        liTY3 = liTY1 * liTY2; 
        liTY4 = liTY1 * liTY3; 
        liTY5 = liTY1 * liTY4;

        /* Compute interpolated value */
        liIV = liVC[ 0]                 + 
               liVC[ 1] * liTY1         + 
               liVC[ 2] * liTY2         + 
               liVC[ 3] * liTY3         + 
               liVC[ 4] * liTY4         + 
               liVC[ 5] * liTY5         +
               liVC[ 6] * liTX1         + 
               liVC[ 7] * liTY1 * liTX1 + 
               liVC[ 8] * liTY2 * liTX1 + 
               liVC[ 9] * liTY3 * liTX1 + 
               liVC[10] * liTY4 * liTX1 + 
               liVC[11] * liTY5 * liTX1 +
               liVC[12] * liTX2         + 
               liVC[13] * liTY1 * liTX2 + 
               liVC[14] * liTY2 * liTX2 + 
               liVC[15] * liTY3 * liTX2 + 
               liVC[16] * liTY4 * liTX2 + 
               liVC[17] * liTY5 * liTX2 +
               liVC[18] * liTX3         + 
               liVC[19] * liTY1 * liTX3 + 
               liVC[20] * liTY2 * liTX3 + 
               liVC[21] * liTY3 * liTX3 + 
               liVC[22] * liTY4 * liTX3 + 
               liVC[23] * liTY5 * liTX3 +
               liVC[24] * liTX4         + 
               liVC[25] * liTY1 * liTX4 + 
               liVC[26] * liTY2 * liTX4 + 
               liVC[27] * liTY3 * liTX4 + 
               liVC[28] * liTY4 * liTX4 + 
               liVC[29] * liTY5 * liTX4 +
               liVC[30] * liTX5         + 
               liVC[31] * liTY1 * liTX5 + 
               liVC[32] * liTY2 * liTX5 + 
               liVC[33] * liTY3 * liTX5 + 
               liVC[34] * liTY4 * liTX5 + 
               liVC[35] * liTY5 * liTX5;

        /* Verify interpolated value */
        liIV = ( liIV < li_Real_s(   0.0 ) ) ? li_Real_s(   0.0 ) : liIV; 
        liIV = ( liIV > li_Real_s( 255.0 ) ) ? li_Real_s( 255.0 ) : liIV;

        /* Return interpolated value */
        return( liIV );

    }

