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

    li_C8_t li_bihepticf(

        li_C8_t * liBytes, 
        li_Size_t liWidth,
        li_Size_t liHeight,
        li_Size_t liLayer, 
        li_Size_t liChannel,
        li_Real_t liX,
        li_Real_t liY

    ) {

        /* Interpolation vectors */
        li_Real_t liVS[64] = { li_Real_s( 0.0 ) };
        li_Real_t liVC[64] = { li_Real_s( 0.0 ) };

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
        li_Real_t liTX6 = li_Real_s( 0.0 );
        li_Real_t liTY6 = li_Real_s( 0.0 );
        li_Real_t liTX7 = li_Real_s( 0.0 );
        li_Real_t liTY7 = li_Real_s( 0.0 );

        /* Interpolation variables */
        li_Size_t liPX = li_Size_s( 0 );
        li_Size_t liPY = li_Size_s( 0 );

        /* Sampling variables */
        li_Size_t liPXm3 = li_Size_s( 0 );
        li_Size_t liPXm2 = li_Size_s( 0 );
        li_Size_t liPXm1 = li_Size_s( 0 );
        li_Size_t liPXp1 = li_Size_s( 0 );
        li_Size_t liPXp2 = li_Size_s( 0 );
        li_Size_t liPXp3 = li_Size_s( 0 );
        li_Size_t liPXp4 = li_Size_s( 0 );
        li_Size_t liPYm3 = li_Size_s( 0 );
        li_Size_t liPYm2 = li_Size_s( 0 );
        li_Size_t liPYm1 = li_Size_s( 0 );
        li_Size_t liPYp1 = li_Size_s( 0 );
        li_Size_t liPYp2 = li_Size_s( 0 );
        li_Size_t liPYp3 = li_Size_s( 0 );
        li_Size_t liPYp4 = li_Size_s( 0 );

        /* Interpolated variables */
        li_Real_t liIV = li_Real_s( 0.0 );

        /* Compute relatliIVe grid parameters */
        liPX = li_Trunc( liX );
        liPY = li_Trunc( liY );

        /* Compute sampling nodes */
        liPXp1 = liPX + li_Size_s( 1 );
        liPYp1 = liPY + li_Size_s( 1 );

        /* Compute sampling nodes */
        liPXm3 = liPX - li_Size_s( 3 ); liPXm3 = ( ( liPXm3 <  li_Size_s( 0 ) ) ? li_Size_s( 0 ) : liPXm3 );
        liPXm2 = liPX - li_Size_s( 2 ); liPXm2 = ( ( liPXm2 <  li_Size_s( 0 ) ) ? li_Size_s( 0 ) : liPXm2 );
        liPXm1 = liPX - li_Size_s( 1 ); liPXm1 = ( ( liPXm1 <  li_Size_s( 0 ) ) ? li_Size_s( 0 ) : liPXm1 );
        liPYm3 = liPY - li_Size_s( 3 ); liPYm3 = ( ( liPYm3 <  li_Size_s( 0 ) ) ? li_Size_s( 0 ) : liPYm3 );
        liPYm2 = liPY - li_Size_s( 2 ); liPYm2 = ( ( liPYm2 <  li_Size_s( 0 ) ) ? li_Size_s( 0 ) : liPYm2 );
        liPYm1 = liPY - li_Size_s( 1 ); liPYm1 = ( ( liPYm1 <  li_Size_s( 0 ) ) ? li_Size_s( 0 ) : liPYm1 );

        /* Compute sampling nodes */
        liPXp2 = liPX + li_Size_s( 2 ); liPXp2 = ( ( liPXp2 >= liWidth  ) ? liWidth  - li_Size_s( 1 ) : liPXp2 );
        liPXp3 = liPX + li_Size_s( 3 ); liPXp3 = ( ( liPXp3 >= liWidth  ) ? liWidth  - li_Size_s( 1 ) : liPXp3 );
        liPXp4 = liPX + li_Size_s( 4 ); liPXp4 = ( ( liPXp4 >= liWidth  ) ? liWidth  - li_Size_s( 1 ) : liPXp4 );
        liPYp2 = liPY + li_Size_s( 2 ); liPYp2 = ( ( liPYp2 >= liHeight ) ? liHeight - li_Size_s( 1 ) : liPYp2 );
        liPYp3 = liPY + li_Size_s( 3 ); liPYp3 = ( ( liPYp3 >= liHeight ) ? liHeight - li_Size_s( 1 ) : liPYp3 );
        liPYp4 = liPY + li_Size_s( 4 ); liPYp4 = ( ( liPYp4 >= liHeight ) ? liHeight - li_Size_s( 1 ) : liPYp4 );

        /* Compute interpolation vector */
        liVS[ 0] = * ( liBytes + liLayer * ( liWidth * liPYm3 + liPXm3 ) + liChannel );
        liVS[ 1] = * ( liBytes + liLayer * ( liWidth * liPYm3 + liPXm2 ) + liChannel );
        liVS[ 2] = * ( liBytes + liLayer * ( liWidth * liPYm3 + liPXm1 ) + liChannel );
        liVS[ 3] = * ( liBytes + liLayer * ( liWidth * liPYm3 + liPX   ) + liChannel );
        liVS[ 4] = * ( liBytes + liLayer * ( liWidth * liPYm3 + liPXp1 ) + liChannel );
        liVS[ 5] = * ( liBytes + liLayer * ( liWidth * liPYm3 + liPXp2 ) + liChannel );
        liVS[ 6] = * ( liBytes + liLayer * ( liWidth * liPYm3 + liPXp3 ) + liChannel );
        liVS[ 7] = * ( liBytes + liLayer * ( liWidth * liPYm3 + liPXp4 ) + liChannel );
        liVS[ 8] = * ( liBytes + liLayer * ( liWidth * liPYm2 + liPXm3 ) + liChannel );
        liVS[ 9] = * ( liBytes + liLayer * ( liWidth * liPYm2 + liPXm2 ) + liChannel );
        liVS[10] = * ( liBytes + liLayer * ( liWidth * liPYm2 + liPXm1 ) + liChannel );
        liVS[11] = * ( liBytes + liLayer * ( liWidth * liPYm2 + liPX   ) + liChannel );
        liVS[12] = * ( liBytes + liLayer * ( liWidth * liPYm2 + liPXp1 ) + liChannel );
        liVS[13] = * ( liBytes + liLayer * ( liWidth * liPYm2 + liPXp2 ) + liChannel );
        liVS[14] = * ( liBytes + liLayer * ( liWidth * liPYm2 + liPXp3 ) + liChannel );
        liVS[15] = * ( liBytes + liLayer * ( liWidth * liPYm2 + liPXp4 ) + liChannel );
        liVS[16] = * ( liBytes + liLayer * ( liWidth * liPYm1 + liPXm3 ) + liChannel );
        liVS[17] = * ( liBytes + liLayer * ( liWidth * liPYm1 + liPXm2 ) + liChannel );
        liVS[18] = * ( liBytes + liLayer * ( liWidth * liPYm1 + liPXm1 ) + liChannel );
        liVS[19] = * ( liBytes + liLayer * ( liWidth * liPYm1 + liPX   ) + liChannel );
        liVS[20] = * ( liBytes + liLayer * ( liWidth * liPYm1 + liPXp1 ) + liChannel );
        liVS[21] = * ( liBytes + liLayer * ( liWidth * liPYm1 + liPXp2 ) + liChannel );
        liVS[22] = * ( liBytes + liLayer * ( liWidth * liPYm1 + liPXp3 ) + liChannel );
        liVS[23] = * ( liBytes + liLayer * ( liWidth * liPYm1 + liPXp4 ) + liChannel );
        liVS[24] = * ( liBytes + liLayer * ( liWidth * liPY   + liPXm3 ) + liChannel );
        liVS[25] = * ( liBytes + liLayer * ( liWidth * liPY   + liPXm2 ) + liChannel );
        liVS[26] = * ( liBytes + liLayer * ( liWidth * liPY   + liPXm1 ) + liChannel );
        liVS[27] = * ( liBytes + liLayer * ( liWidth * liPY   + liPX   ) + liChannel );
        liVS[28] = * ( liBytes + liLayer * ( liWidth * liPY   + liPXp1 ) + liChannel );
        liVS[29] = * ( liBytes + liLayer * ( liWidth * liPY   + liPXp2 ) + liChannel );
        liVS[30] = * ( liBytes + liLayer * ( liWidth * liPY   + liPXp3 ) + liChannel );
        liVS[31] = * ( liBytes + liLayer * ( liWidth * liPY   + liPXp4 ) + liChannel );
        liVS[32] = * ( liBytes + liLayer * ( liWidth * liPYp1 + liPXm3 ) + liChannel );
        liVS[33] = * ( liBytes + liLayer * ( liWidth * liPYp1 + liPXm2 ) + liChannel );
        liVS[34] = * ( liBytes + liLayer * ( liWidth * liPYp1 + liPXm1 ) + liChannel );
        liVS[35] = * ( liBytes + liLayer * ( liWidth * liPYp1 + liPX   ) + liChannel );
        liVS[36] = * ( liBytes + liLayer * ( liWidth * liPYp1 + liPXp1 ) + liChannel );
        liVS[37] = * ( liBytes + liLayer * ( liWidth * liPYp1 + liPXp2 ) + liChannel );
        liVS[38] = * ( liBytes + liLayer * ( liWidth * liPYp1 + liPXp3 ) + liChannel );
        liVS[39] = * ( liBytes + liLayer * ( liWidth * liPYp1 + liPXp4 ) + liChannel );
        liVS[40] = * ( liBytes + liLayer * ( liWidth * liPYp2 + liPXm3 ) + liChannel );
        liVS[41] = * ( liBytes + liLayer * ( liWidth * liPYp2 + liPXm2 ) + liChannel );
        liVS[42] = * ( liBytes + liLayer * ( liWidth * liPYp2 + liPXm1 ) + liChannel );
        liVS[43] = * ( liBytes + liLayer * ( liWidth * liPYp2 + liPX   ) + liChannel );
        liVS[44] = * ( liBytes + liLayer * ( liWidth * liPYp2 + liPXp1 ) + liChannel );
        liVS[45] = * ( liBytes + liLayer * ( liWidth * liPYp2 + liPXp2 ) + liChannel );
        liVS[46] = * ( liBytes + liLayer * ( liWidth * liPYp2 + liPXp3 ) + liChannel );
        liVS[47] = * ( liBytes + liLayer * ( liWidth * liPYp2 + liPXp4 ) + liChannel );
        liVS[48] = * ( liBytes + liLayer * ( liWidth * liPYp3 + liPXm3 ) + liChannel );
        liVS[49] = * ( liBytes + liLayer * ( liWidth * liPYp3 + liPXm2 ) + liChannel );
        liVS[50] = * ( liBytes + liLayer * ( liWidth * liPYp3 + liPXm1 ) + liChannel );
        liVS[51] = * ( liBytes + liLayer * ( liWidth * liPYp3 + liPX   ) + liChannel );
        liVS[52] = * ( liBytes + liLayer * ( liWidth * liPYp3 + liPXp1 ) + liChannel );
        liVS[53] = * ( liBytes + liLayer * ( liWidth * liPYp3 + liPXp2 ) + liChannel );
        liVS[54] = * ( liBytes + liLayer * ( liWidth * liPYp3 + liPXp3 ) + liChannel );
        liVS[55] = * ( liBytes + liLayer * ( liWidth * liPYp3 + liPXp4 ) + liChannel );
        liVS[56] = * ( liBytes + liLayer * ( liWidth * liPYp4 + liPXm3 ) + liChannel );
        liVS[57] = * ( liBytes + liLayer * ( liWidth * liPYp4 + liPXm2 ) + liChannel );
        liVS[58] = * ( liBytes + liLayer * ( liWidth * liPYp4 + liPXm1 ) + liChannel );
        liVS[59] = * ( liBytes + liLayer * ( liWidth * liPYp4 + liPX   ) + liChannel );
        liVS[60] = * ( liBytes + liLayer * ( liWidth * liPYp4 + liPXp1 ) + liChannel );
        liVS[61] = * ( liBytes + liLayer * ( liWidth * liPYp4 + liPXp2 ) + liChannel );
        liVS[62] = * ( liBytes + liLayer * ( liWidth * liPYp4 + liPXp3 ) + liChannel );
        liVS[63] = * ( liBytes + liLayer * ( liWidth * liPYp4 + liPXp4 ) + liChannel );

        /* Compute interpolation matrix product */
        liVC[ 0] = li_Real_s(       1.000000000000000000000000 ) * liVS[ 0];
        liVC[ 1] = li_Real_s(      -2.592857142857315277950647 ) * liVS[ 0] +
                   li_Real_s(       7.000000000001223021683927 ) * liVS[ 8] +
                   li_Real_s(     -10.500000000003735678433259 ) * liVS[16] +
                   li_Real_s(      11.666666666672959706829715 ) * liVS[24] +
                   li_Real_s(      -8.750000000006311395850389 ) * liVS[32] +
                   li_Real_s(       4.200000000003775824097829 ) * liVS[40] +
                   li_Real_s(      -1.166666666667915297495028 ) * liVS[48] +
                   li_Real_s(       0.142857142857318653028642 ) * liVS[56];
        liVC[ 2] = li_Real_s(       2.605555555556023250574071 ) * liVS[ 0] +
                   li_Real_s(     -11.150000000003235101075916 ) * liVS[ 8] +
                   li_Real_s(      21.975000000009650591437094 ) * liVS[16] +
                   li_Real_s(     -26.361111111127101480633428 ) * liVS[24] +
                   li_Real_s(      20.500000000015841550293771 ) * liVS[32] +
                   li_Real_s(     -10.050000000009386980082127 ) * liVS[40] +
                   li_Real_s(       2.830555555558634139856622 ) * liVS[48] +
                   li_Real_s(      -0.350000000000429523083767 ) * liVS[56];
        liVC[ 3] = li_Real_s(      -1.343055555555963564984268 ) * liVS[ 0] +
                   li_Real_s(       7.088888888891666795188939 ) * liVS[ 8] +
                   li_Real_s(     -16.370833333341497706214795 ) * liVS[16] +
                   li_Real_s(      21.611111111124493788793188 ) * liVS[24] +
                   li_Real_s(     -17.673611111124262862404066 ) * liVS[32] +
                   li_Real_s(       8.933333333341073156930179 ) * liVS[40] +
                   li_Real_s(      -2.568055555558079738887045 ) * liVS[48] +
                   li_Real_s(       0.322222222222574572469966 ) * liVS[56];
        liVC[ 4] = li_Real_s(       0.388888888889063366605114 ) * liVS[ 0] +
                   li_Real_s(      -2.312500000001146638339833 ) * liVS[ 8] +
                   li_Real_s(       5.916666666669997631800015 ) * liVS[16] +
                   li_Real_s(      -8.465277777783192902916198 ) * liVS[24] +
                   li_Real_s(       7.333333333338621251584755 ) * liVS[32] +
                   li_Real_s(      -3.854166666669763596786424 ) * liVS[40] +
                   li_Real_s(       1.138888888889894479561349 ) * liVS[48] +
                   li_Real_s(      -0.145833333333472481285753 ) * liVS[56];
        liVC[ 5] = li_Real_s(      -0.063888888888927075626611 ) * liVS[ 0] +
                   li_Real_s(       0.409722222222470067176658 ) * liVS[ 8] +
                   li_Real_s(      -1.125000000000712985226414 ) * liVS[16] +
                   li_Real_s(       1.715277777778928314234008 ) * liVS[24] +
                   li_Real_s(      -1.569444444445561748224804 ) * liVS[32] +
                   li_Real_s(       0.862500000000651301235166 ) * liVS[40] +
                   li_Real_s(      -0.263888888889099559875717 ) * liVS[48] +
                   li_Real_s(       0.034722222222251408751958 ) * liVS[56];
        liVC[ 6] = li_Real_s(       0.005555555555559685276812 ) * liVS[ 0] +
                   li_Real_s(      -0.037500000000026852131629 ) * liVS[ 8] +
                   li_Real_s(       0.108333333333409928544988 ) * liVS[16] +
                   li_Real_s(      -0.173611111111233978876456 ) * liVS[24] +
                   li_Real_s(       0.166666666666785423522867 ) * liVS[32] +
                   li_Real_s(      -0.095833333333402298537251 ) * liVS[40] +
                   li_Real_s(       0.030555555555577783299892 ) * liVS[48] +
                   li_Real_s(      -0.004166666666669732732586 ) * liVS[56];
        liVC[ 7] = li_Real_s(      -0.000198412698412875516951 ) * liVS[ 0] +
                   li_Real_s(       0.001388888888890034076229 ) * liVS[ 8] +
                   li_Real_s(      -0.004166666666669913143828 ) * liVS[16] +
                   li_Real_s(       0.006944444444449626545335 ) * liVS[24] +
                   li_Real_s(      -0.006944444444449433991029 ) * liVS[32] +
                   li_Real_s(       0.004166666666669554056068 ) * liVS[40] +
                   li_Real_s(      -0.001388888888889816802114 ) * liVS[48] +
                   li_Real_s(       0.000198412698412826077332 ) * liVS[56];
        liVC[ 8] = li_Real_s(      -2.592857143098429073546640 ) * liVS[ 0] +
                   li_Real_s(       7.000000001583480901956591 ) * liVS[ 1] +
                   li_Real_s(     -10.500000004488029503590951 ) * liVS[ 2] +
                   li_Real_s(      11.666666673745227811309633 ) * liVS[ 3] +
                   li_Real_s(      -8.750000006691607268294320 ) * liVS[ 4] +
                   li_Real_s(       4.200000003790139047055163 ) * liVS[ 5] +
                   li_Real_s(      -1.166666667858553108061415 ) * liVS[ 6] +
                   li_Real_s(       0.142857143017770776838304 ) * liVS[ 7] +
                   li_Real_s(       0.000000001540405075126256 ) * liVS[ 8] +
                   li_Real_s(      -0.000000010136165283804684 ) * liVS[ 9] +
                   li_Real_s(       0.000000028755684795670651 ) * liVS[10] +
                   li_Real_s(      -0.000000045346676490054424 ) * liVS[11] +
                   li_Real_s(       0.000000042835230190816962 ) * liVS[12] +
                   li_Real_s(      -0.000000024237404063520649 ) * liVS[13] +
                   li_Real_s(       0.000000007614139591219153 ) * liVS[14] +
                   li_Real_s(      -0.000000001025213815453246 ) * liVS[15] +
                   li_Real_s(      -0.000000004186690509776294 ) * liVS[16] +
                   li_Real_s(       0.000000027631880062657426 ) * liVS[17] +
                   li_Real_s(      -0.000000078492344722262946 ) * liVS[18] +
                   li_Real_s(       0.000000123799260595672108 ) * liVS[19] +
                   li_Real_s(      -0.000000116880915515000991 ) * liVS[20] +
                   li_Real_s(       0.000000066077549464921399 ) * liVS[21] +
                   li_Real_s(      -0.000000020738669988760520 ) * liVS[22] +
                   li_Real_s(       0.000000002789930612549794 ) * liVS[23] +
                   li_Real_s(       0.000000006264439020481958 ) * liVS[24] +
                   li_Real_s(      -0.000000041485519063464229 ) * liVS[25] +
                   li_Real_s(       0.000000118053501280153854 ) * liVS[26] +
                   li_Real_s(      -0.000000186303181648561394 ) * liVS[27] +
                   li_Real_s(       0.000000175858372337865625 ) * liVS[28] +
                   li_Real_s(      -0.000000099359597316023034 ) * liVS[29] +
                   li_Real_s(       0.000000031160779741088901 ) * liVS[30] +
                   li_Real_s(      -0.000000004188794351541618 ) * liVS[31] +
                   li_Real_s(      -0.000000005579556117742455 ) * liVS[32] +
                   li_Real_s(       0.000000037090397610960538 ) * liVS[33] +
                   li_Real_s(      -0.000000105781213827561252 ) * liVS[34] +
                   li_Real_s(       0.000000167104460118753317 ) * liVS[35] +
                   li_Real_s(      -0.000000157763520204439358 ) * liVS[36] +
                   li_Real_s(       0.000000089107847277009771 ) * liVS[37] +
                   li_Real_s(      -0.000000027930730990698894 ) * liVS[38] +
                   li_Real_s(       0.000000003752316133718326 ) * liVS[39] +
                   li_Real_s(       0.000000002966981750542554 ) * liVS[40] +
                   li_Real_s(      -0.000000019804072490595199 ) * liVS[41] +
                   li_Real_s(       0.000000056627451408063459 ) * liVS[42] +
                   li_Real_s(      -0.000000089578011431070257 ) * liVS[43] +
                   li_Real_s(       0.000000084612122765756188 ) * liVS[44] +
                   li_Real_s(      -0.000000047787756143580296 ) * liVS[45] +
                   li_Real_s(       0.000000014974042047580831 ) * liVS[46] +
                   li_Real_s(      -0.000000002010757906697286 ) * liVS[47] +
                   li_Real_s(      -0.000000000875486708153441 ) * liVS[48] +
                   li_Real_s(       0.000000005868271200811676 ) * liVS[49] +
                   li_Real_s(      -0.000000016826386810402470 ) * liVS[50] +
                   li_Real_s(       0.000000026659626149677877 ) * liVS[51] +
                   li_Real_s(      -0.000000025199213124527652 ) * liVS[52] +
                   li_Real_s(       0.000000014233864279363162 ) * liVS[53] +
                   li_Real_s(      -0.000000004459275001327134 ) * liVS[54] +
                   li_Real_s(       0.000000000598600014557984 ) * liVS[55] +
                   li_Real_s(       0.000000000111060541230283 ) * liVS[56] +
                   li_Real_s(      -0.000000000747433628892725 ) * liVS[57] +
                   li_Real_s(       0.000000002149041305809822 ) * liVS[58] +
                   li_Real_s(      -0.000000003410454989932734 ) * liVS[59] +
                   li_Real_s(       0.000000003226093757437956 ) * liVS[60] +
                   li_Real_s(      -0.000000001822630869507494 ) * liVS[61] +
                   li_Real_s(       0.000000000570939797403020 ) * liVS[62] +
                   li_Real_s(      -0.000000000076615913548116 ) * liVS[63];
        liVC[ 9] = li_Real_s(       6.722908163561157834919868 ) * liVS[ 0] +
                   li_Real_s(     -18.150000001730816734379914 ) * liVS[ 1] +
                   li_Real_s(      27.225000004749126958358829 ) * liVS[ 2] +
                   li_Real_s(     -30.250000007578194072266342 ) * liVS[ 3] +
                   li_Real_s(      22.687500007316316441574600 ) * liVS[ 4] +
                   li_Real_s(     -10.890000004182780912742601 ) * liVS[ 5] +
                   li_Real_s(       3.025000001298232543689437 ) * liVS[ 6] +
                   li_Real_s(      -0.370408163433040726886247 ) * liVS[ 7] +
                   li_Real_s(     -18.150000003363722100857558 ) * liVS[ 8] +
                   li_Real_s(      49.000000021454027887557459 ) * liVS[ 9] +
                   li_Real_s(     -73.500000061087149560989928 ) * liVS[10] +
                   li_Real_s(      81.666666764825237123659463 ) * liVS[11] +
                   li_Real_s(     -61.250000094498105340790062 ) * liVS[12] +
                   li_Real_s(      29.400000054013279537912240 ) * liVS[13] +
                   li_Real_s(      -8.166666683579496321954139 ) * liVS[14] +
                   li_Real_s(       1.000000002235935880889883 ) * liVS[15] +
                   li_Real_s(      27.225000012390047743338073 ) * liVS[16] +
                   li_Real_s(     -73.500000081052306200035673 ) * liVS[17] +
                   li_Real_s(     110.250000232936613997480890 ) * liVS[18] +
                   li_Real_s(    -122.500000374620412912918255 ) * liVS[19] +
                   li_Real_s(      91.875000360039166480419226 ) * liVS[20] +
                   li_Real_s(     -44.100000205602832181739359 ) * liVS[21] +
                   li_Real_s(      12.250000064478090422426249 ) * liVS[22] +
                   li_Real_s(      -1.500000008568385112539545 ) * liVS[23] +
                   li_Real_s(     -30.250000022284737610789307 ) * liVS[24] +
                   li_Real_s(      81.666666814221613890367735 ) * liVS[25] +
                   li_Real_s(    -122.500000425965041017661861 ) * liVS[26] +
                   li_Real_s(     136.111111796570895648983424 ) * liVS[27] +
                   li_Real_s(    -102.083333991649752192643064 ) * liVS[28] +
                   li_Real_s(      49.000000375732803092887480 ) * liVS[29] +
                   li_Real_s(     -13.611111228984633214622590 ) * liVS[30] +
                   li_Real_s(       1.666666682358838302846493 ) * liVS[31] +
                   li_Real_s(      22.687500022489771822620241 ) * liVS[32] +
                   li_Real_s(     -61.250000150042012592166429 ) * liVS[33] +
                   li_Real_s(      91.875000434443208519041946 ) * liVS[34] +
                   li_Real_s(    -102.083334032869771590412711 ) * liVS[35] +
                   li_Real_s(      76.562500671639725169370649 ) * liVS[36] +
                   li_Real_s(     -36.750000383182353402844456 ) * liVS[37] +
                   li_Real_s(      10.208333453531633239208531 ) * liVS[38] +
                   li_Real_s(      -1.250000016010197612104093 ) * liVS[39] +
                   li_Real_s(     -10.890000013115809451846872 ) * liVS[40] +
                   li_Real_s(      29.400000087990868280485302 ) * liVS[41] +
                   li_Real_s(     -44.100000255383918101870222 ) * liVS[42] +
                   li_Real_s(      49.000000411495221896984731 ) * liVS[43] +
                   li_Real_s(     -36.750000395049994494911516 ) * liVS[44] +
                   li_Real_s(      17.640000225305339398573778 ) * liVS[45] +
                   li_Real_s(      -4.900000070651838512958420 ) * liVS[46] +
                   li_Real_s(       0.600000009410102563833789 ) * liVS[47] +
                   li_Real_s(       3.025000004153689303620922 ) * liVS[48] +
                   li_Real_s(      -8.166666694658744063417544 ) * liVS[49] +
                   li_Real_s(      12.250000081414796682111046 ) * liVS[50] +
                   li_Real_s(     -13.611111242391643827431835 ) * liVS[51] +
                   li_Real_s(      10.208333459374557605769951 ) * liVS[52] +
                   li_Real_s(      -4.900000071863917838754787 ) * liVS[53] +
                   li_Real_s(       1.361111133637109560368117 ) * liVS[54] +
                   li_Real_s(      -0.166666669665840316838512 ) * liVS[55] +
                   li_Real_s(      -0.370408163820314939584932 ) * liVS[56] +
                   li_Real_s(       1.000000003754809574374462 ) * liVS[57] +
                   li_Real_s(      -1.500000010941763051164344 ) * liVS[58] +
                   li_Real_s(       1.666666684324450642407101 ) * liVS[59] +
                   li_Real_s(      -1.250000016955883808122962 ) * liVS[60] +
                   li_Real_s(       0.600000009665480171783258 ) * liVS[61] +
                   li_Real_s(      -0.166666669694976121718355 ) * liVS[62] +
                   li_Real_s(       0.020408163668268031187836 ) * liVS[63];
        liVC[10] = li_Real_s(      -6.755833332806446378526744 ) * liVS[ 0] +
                   li_Real_s(      18.238888885188437427586905 ) * liVS[ 1] +
                   li_Real_s(     -27.358333323148006144265310 ) * liVS[ 2] +
                   li_Real_s(      30.398148133213396704377374 ) * liVS[ 3] +
                   li_Real_s(     -22.798611098079781811520661 ) * liVS[ 4] +
                   li_Real_s(      10.943333326379345749046479 ) * liVS[ 5] +
                   li_Real_s(      -3.039814812673137467413653 ) * liVS[ 6] +
                   li_Real_s(       0.372222221926246987777631 ) * liVS[ 7] +
                   li_Real_s(      28.910357142740366498401272 ) * liVS[ 8] +
                   li_Real_s(     -78.049999999230266212180140 ) * liVS[ 9] +
                   li_Real_s(     117.075000002556905087658379 ) * liVS[10] +
                   li_Real_s(    -130.083333347429032755826483 ) * liVS[11] +
                   li_Real_s(      97.562500022021538370609051 ) * liVS[12] +
                   li_Real_s(     -46.830000015818917802334909 ) * liVS[13] +
                   li_Real_s(      13.008333338699458181508817 ) * liVS[14] +
                   li_Real_s(      -1.592857143540044262408628 ) * liVS[15] +
                   li_Real_s(     -56.978035720978553513305087 ) * liVS[16] +
                   li_Real_s(     153.825000046907774731153040 ) * liVS[17] +
                   li_Real_s(    -230.737500151161185613091220 ) * liVS[18] +
                   li_Real_s(     256.375000271107126081915339 ) * liVS[19] +
                   li_Real_s(    -192.281250283115554111645906 ) * liVS[20] +
                   li_Real_s(      92.295000170460170352271234 ) * liVS[21] +
                   li_Real_s(     -25.637500054777341063072527 ) * liVS[22] +
                   li_Real_s(       3.139285721557541819493053 ) * liVS[23] +
                   li_Real_s(      68.350595255907222735913820 ) * liVS[24] +
                   li_Real_s(    -184.527777902619334327027900 ) * liVS[25] +
                   li_Real_s(     276.791667054066806485934649 ) * liVS[26] +
                   li_Real_s(    -307.546296962698647803335916 ) * liVS[27] +
                   li_Real_s(     230.659722896195916064243647 ) * liVS[28] +
                   li_Real_s(    -110.716667064636510531272506 ) * liVS[29] +
                   li_Real_s(      30.754629756582986033208726 ) * liVS[30] +
                   li_Real_s(      -3.765873032798785491337412 ) * liVS[31] +
                   li_Real_s(     -53.153571449653043146099662 ) * liVS[32] +
                   li_Real_s(     143.500000148003522326689563 ) * liVS[33] +
                   li_Real_s(    -215.250000454173488151354832 ) * liVS[34] +
                   li_Real_s(     239.166667437337110868611489 ) * liVS[35] +
                   li_Real_s(    -179.375000770667469396357774 ) * liVS[36] +
                   li_Real_s(      86.100000451698946335454821 ) * liVS[37] +
                   li_Real_s(     -23.916666810277828147945911 ) * liVS[38] +
                   li_Real_s(       2.928571447732210231151839 ) * liVS[39] +
                   li_Real_s(      26.058214299061674523727561 ) * liVS[40] +
                   li_Real_s(     -70.350000093943265255802544 ) * liVS[41] +
                   li_Real_s(     105.525000287023786427198502 ) * liVS[42] +
                   li_Real_s(    -117.250000484044647919290583 ) * liVS[43] +
                   li_Real_s(      87.937500481389363926609803 ) * liVS[44] +
                   li_Real_s(     -42.210000281035512159633072 ) * liVS[45] +
                   li_Real_s(      11.725000089153226667804120 ) * liVS[46] +
                   li_Real_s(      -1.435714297604568034927297 ) * liVS[47] +
                   li_Real_s(      -7.339226194896596666694677 ) * liVS[48] +
                   li_Real_s(      19.813888920104236746055903 ) * liVS[49] +
                   li_Real_s(     -29.720833428579311430439702 ) * liVS[50] +
                   li_Real_s(      33.023148308314603127655573 ) * liVS[51] +
                   li_Real_s(     -24.767361269939499379688641 ) * liVS[52] +
                   li_Real_s(      11.888333425839512358379579 ) * liVS[53] +
                   li_Real_s(      -3.302314844112274272447394 ) * liVS[54] +
                   li_Real_s(       0.404365083269311753610964 ) * liVS[55] +
                   li_Real_s(       0.907500000603647549723973 ) * liVS[56] +
                   li_Real_s(      -2.450000004276803089453551 ) * liVS[57] +
                   li_Real_s(       3.675000013056727965476966 ) * liVS[58] +
                   li_Real_s(      -4.083333355270326592290075 ) * liVS[59] +
                   li_Real_s(       3.062500021725167442809834 ) * liVS[60] +
                   li_Real_s(      -1.470000012636015762623742 ) * liVS[61] +
                   li_Real_s(       0.408333337330143208987465 ) * liVS[62] +
                   li_Real_s(      -0.050000000532236299477518 ) * liVS[63];
        liVC[11] = li_Real_s(       3.482351189884582254308043 ) * liVS[ 0] +
                   li_Real_s(      -9.401388884821457736507000 ) * liVS[ 1] +
                   li_Real_s(      14.102083322077589855325641 ) * liVS[ 2] +
                   li_Real_s(     -15.668981464666217107151169 ) * liVS[ 3] +
                   li_Real_s(      11.751736096081238258648227 ) * liVS[ 4] +
                   li_Real_s(      -5.640833325132277309421625 ) * liVS[ 5] +
                   li_Real_s(       1.566898145586483437341485 ) * liVS[ 6] +
                   li_Real_s(      -0.191865079010135275439097 ) * liVS[ 7] +
                   li_Real_s(     -18.380476189374412854249385 ) * liVS[ 8] +
                   li_Real_s(      49.622222215124111244222149 ) * liVS[ 9] +
                   li_Real_s(     -74.433333317399430484329059 ) * liVS[10] +
                   li_Real_s(      82.703703687112579245876987 ) * liVS[11] +
                   li_Real_s(     -62.027777769257262718838319 ) * liVS[12] +
                   li_Real_s(      29.773333331104698373792417 ) * liVS[13] +
                   li_Real_s(      -8.270370369891647754911901 ) * liVS[14] +
                   li_Real_s(       1.012698412581357843009755 ) * liVS[15] +
                   li_Real_s(      42.447232145623857491045783 ) * liVS[16] +
                   li_Real_s(    -114.595833354374121881846804 ) * liVS[17] +
                   li_Real_s(     171.893750075601133175950963 ) * liVS[18] +
                   li_Real_s(    -190.993055702750893942720722 ) * liVS[19] +
                   li_Real_s(     143.244791828234980357592576 ) * liVS[20] +
                   li_Real_s(     -68.757500099563941375890863 ) * liVS[21] +
                   li_Real_s(      19.099305587598912836710952 ) * liVS[22] +
                   li_Real_s(      -2.338690480369933766269241 ) * liVS[23] +
                   li_Real_s(     -56.034523820060641696727544 ) * liVS[24] +
                   li_Real_s(     151.277777854103646859584842 ) * liVS[25] +
                   li_Real_s(    -226.916666912157779734116048 ) * liVS[26] +
                   li_Real_s(     252.129630064173056780418847 ) * liVS[27] +
                   li_Real_s(    -189.097222669932648386748042 ) * liVS[28] +
                   li_Real_s(      90.766666933306481723775505 ) * liVS[29] +
                   li_Real_s(     -25.212963047967015484118747 ) * liVS[30] +
                   li_Real_s(       3.087301598534923474659308 ) * liVS[31] +
                   li_Real_s(      45.825148823163260658475338 ) * liVS[32] +
                   li_Real_s(    -123.715277875749634972635249 ) * liVS[33] +
                   li_Real_s(     185.572916974194754402560648 ) * liVS[34] +
                   li_Real_s(    -206.192130160933146498791757 ) * liVS[35] +
                   li_Real_s(     154.644097759783278434042586 ) * liVS[36] +
                   li_Real_s(     -74.229166983363001008910942 ) * liVS[37] +
                   li_Real_s(      20.619213063533397445326045 ) * liVS[38] +
                   li_Real_s(      -2.524801600628919118207705 ) * liVS[39] +
                   li_Real_s(     -23.162857151830706925466075 ) * liVS[40] +
                   li_Real_s(      62.533333397683456666982238 ) * liVS[41] +
                   li_Real_s(     -93.800000200089243662660010 ) * liVS[42] +
                   li_Real_s(     104.222222564351710616392666 ) * liVS[43] +
                   li_Real_s(     -78.166667009968378465600836 ) * liVS[44] +
                   li_Real_s(      37.520000201153322905156529 ) * liVS[45] +
                   li_Real_s(     -10.422222285936012298179776 ) * liVS[46] +
                   li_Real_s(       1.276190484635840505234228 ) * liVS[47] +
                   li_Real_s(       6.658601193502974524562887 ) * liVS[48] +
                   li_Real_s(     -17.976388910624677919258829 ) * liVS[49] +
                   li_Real_s(      26.964583400671187973784981 ) * liVS[50] +
                   li_Real_s(     -29.960648262733428737192298 ) * liVS[51] +
                   li_Real_s(      22.470486225600524221590604 ) * liVS[52] +
                   li_Real_s(     -10.785833400206382037822550 ) * liVS[53] +
                   li_Real_s(       2.996064835955138505596551 ) * liVS[54] +
                   li_Real_s(      -0.366865082165347189402382 ) * liVS[55] +
                   li_Real_s(      -0.835476190892279646504903 ) * liVS[56] +
                   li_Real_s(       2.255555558553805184374141 ) * liVS[57] +
                   li_Real_s(      -3.383333342618197292495097 ) * liVS[58] +
                   li_Real_s(       3.759259275031425318047695 ) * liVS[59] +
                   li_Real_s(      -2.819444460172398692066054 ) * liVS[60] +
                   li_Real_s(       1.353333342503419522984132 ) * liVS[61] +
                   li_Real_s(      -0.375925928820343813185900 ) * liVS[62] +
                   li_Real_s(       0.046031746414648466725339 ) * liVS[63];
        liVC[12] = li_Real_s(      -1.008333333254626040798030 ) * liVS[ 0] +
                   li_Real_s(       2.722222221631007244013745 ) * liVS[ 1] +
                   li_Real_s(      -4.083333331736659133071043 ) * liVS[ 2] +
                   li_Real_s(       4.537037034831897130970901 ) * liVS[ 3] +
                   li_Real_s(      -3.402777775979600960454263 ) * liVS[ 4] +
                   li_Real_s(       1.633333332405122462205327 ) * liVS[ 5] +
                   li_Real_s(      -0.453703703407383862611368 ) * liVS[ 6] +
                   li_Real_s(       0.055555555510220955284240 ) * liVS[ 7] +
                   li_Real_s(       5.995982143324988911103901 ) * liVS[ 8] +
                   li_Real_s(     -16.187500003146009674992456 ) * liVS[ 9] +
                   li_Real_s(      24.281250010472263056726661 ) * liVS[10] +
                   li_Real_s(     -26.979166686312254341828520 ) * liVS[11] +
                   li_Real_s(      20.234375021138504280315829 ) * liVS[12] +
                   li_Real_s(      -9.712500012829925921664653 ) * liVS[13] +
                   li_Real_s(       2.697916670726488774789686 ) * liVS[14] +
                   li_Real_s(      -0.330357143374058637164126 ) * liVS[15] +
                   li_Real_s(     -15.341071431977461259066331 ) * liVS[16] +
                   li_Real_s(      41.416666690272293749330856 ) * liVS[17] +
                   li_Real_s(     -62.125000073094547303753643 ) * liVS[18] +
                   li_Real_s(      69.027777903651283963881724 ) * liVS[19] +
                   li_Real_s(     -51.770833460629923195028823 ) * liVS[20] +
                   li_Real_s(      24.850000074890374435199192 ) * liVS[21] +
                   li_Real_s(      -6.902777801462150364386616 ) * liVS[22] +
                   li_Real_s(       0.845238098350137079250999 ) * liVS[23] +
                   li_Real_s(      21.949255959896611045678583 ) * liVS[24] +
                   li_Real_s(     -59.256944496806610800376802 ) * liVS[25] +
                   li_Real_s(      88.885416826525414535353775 ) * liVS[26] +
                   li_Real_s(     -98.761574344381870105280541 ) * liVS[27] +
                   li_Real_s(      74.071180824959114374905766 ) * liVS[28] +
                   li_Real_s(     -35.554166823868499136551691 ) * liVS[29] +
                   li_Real_s(       9.876157457068213574302717 ) * liVS[30] +
                   li_Real_s(      -1.209325403392290887438776 ) * liVS[31] +
                   li_Real_s(     -19.014285722628301300574094 ) * liVS[32] +
                   li_Real_s(      51.333333391613734875136288 ) * liVS[33] +
                   li_Real_s(     -77.000000176932431372733845 ) * liVS[34] +
                   li_Real_s(      85.555555852453068155227811 ) * liVS[35] +
                   li_Real_s(     -64.166666960659327401117480 ) * liVS[36] +
                   li_Real_s(      30.800000170851188130427545 ) * liVS[37] +
                   li_Real_s(      -8.555555609455982590816348 ) * liVS[38] +
                   li_Real_s(       1.047619054758065715304838 ) * liVS[39] +
                   li_Real_s(       9.993303576548299815840437 ) * liVS[40] +
                   li_Real_s(     -26.979166702519556508832466 ) * liVS[41] +
                   li_Real_s(      40.468750108578930735347967 ) * liVS[42] +
                   li_Real_s(     -44.965277959269087659777142 ) * liVS[43] +
                   li_Real_s(      33.723958512406824183926801 ) * liVS[44] +
                   li_Real_s(     -16.187500103796576667036788 ) * liVS[45] +
                   li_Real_s(       4.496527810477026321223093 ) * liVS[46] +
                   li_Real_s(      -0.550595242425866882030050 ) * liVS[47] +
                   li_Real_s(      -2.952976192140251754381097 ) * liVS[48] +
                   li_Real_s(       7.972222233910155608782588 ) * liVS[49] +
                   li_Real_s(     -11.958333368711308253296011 ) * liVS[50] +
                   li_Real_s(      13.287037096066008245998091 ) * liVS[51] +
                   li_Real_s(      -9.965277835906249492836650 ) * liVS[52] +
                   li_Real_s(       4.783333366968562927468156 ) * liVS[53] +
                   li_Real_s(      -1.328703714285730086430704 ) * liVS[54] +
                   li_Real_s(       0.162698414098812804695626 ) * liVS[55] +
                   li_Real_s(       0.378125000224244445234945 ) * liVS[56] +
                   li_Real_s(      -1.020833334914323486941612 ) * liVS[57] +
                   li_Real_s(       1.531250004789606933286450 ) * liVS[58] +
                   li_Real_s(      -1.701388896877656264905454 ) * liVS[59] +
                   li_Real_s(       1.276041674526961600122377 ) * liVS[60] +
                   li_Real_s(      -0.612500004543345188956494 ) * liVS[61] +
                   li_Real_s(       0.170138890316543722747156 ) * liVS[62] +
                   li_Real_s(      -0.020833333522048746999644 ) * liVS[63];
        liVC[13] = li_Real_s(       0.165654761916623982642705 ) * liVS[ 0] +
                   li_Real_s(      -0.447222222278342407264518 ) * liVS[ 1] +
                   li_Real_s(       0.670833333488994654203452 ) * liVS[ 2] +
                   li_Real_s(      -0.745370370651496472191866 ) * liVS[ 3] +
                   li_Real_s(       0.559027778078156067920190 ) * liVS[ 4] +
                   li_Real_s(      -0.268333333509552574014378 ) * liVS[ 5] +
                   li_Real_s(       0.074537037087859248085664 ) * liVS[ 6] +
                   li_Real_s(      -0.009126984132252324855017 ) * liVS[ 7] +
                   li_Real_s(      -1.062351190757836238987011 ) * liVS[ 8] +
                   li_Real_s(       2.868055557386362419691750 ) * liVS[ 9] +
                   li_Real_s(      -4.302083338736898632248540 ) * liVS[10] +
                   li_Real_s(       4.780092601602168755903222 ) * liVS[11] +
                   li_Real_s(      -3.585069453356973134816599 ) * liVS[12] +
                   li_Real_s(       1.720833338489091479672766 ) * liVS[13] +
                   li_Real_s(      -0.478009260865195573231290 ) * liVS[14] +
                   li_Real_s(       0.058531746239280479926492 ) * liVS[15] +
                   li_Real_s(       2.916964286935748162932214 ) * liVS[16] +
                   li_Real_s(      -7.875000008196489709177968 ) * liVS[17] +
                   li_Real_s(      11.812500024249278141041941 ) * liVS[18] +
                   li_Real_s(     -13.125000040042909432713714 ) * liVS[19] +
                   li_Real_s(       9.843750039238894800064372 ) * liVS[20] +
                   li_Real_s(      -4.725000022630218232677635 ) * liVS[21] +
                   li_Real_s(       1.312500007092642917427838 ) * liVS[22] +
                   li_Real_s(      -0.160714286646944870540210 ) * liVS[23] +
                   li_Real_s(      -4.447470240451119494196064 ) * liVS[24] +
                   li_Real_s(      12.006944460421369313962714 ) * liVS[25] +
                   li_Real_s(     -18.010416713938948163331588 ) * liVS[26] +
                   li_Real_s(      20.011574151824252254527892 ) * liVS[27] +
                   li_Real_s(     -15.008680631445722397643294 ) * liVS[28] +
                   li_Real_s(       7.204166710356727065800442 ) * liVS[29] +
                   li_Real_s(      -2.001157421118117341052312 ) * liVS[30] +
                   li_Real_s(       0.245039684351555570041015 ) * liVS[31] +
                   li_Real_s(       4.069345240565088062112409 ) * liVS[32] +
                   li_Real_s(     -10.986111127946850629655273 ) * liVS[33] +
                   li_Real_s(      16.479166716468043318855052 ) * liVS[34] +
                   li_Real_s(     -18.310185266906611190051990 ) * liVS[35] +
                   li_Real_s(      13.732638968466094553377843 ) * liVS[36] +
                   li_Real_s(      -6.591666712411151429762413 ) * liVS[37] +
                   li_Real_s(       1.831018532871099724346209 ) * liVS[38] +
                   li_Real_s(      -0.224206351105712187177232 ) * liVS[39] +
                   li_Real_s(      -2.236339287187095514752855 ) * liVS[40] +
                   li_Real_s(       6.037500010072765377344695 ) * liVS[41] +
                   li_Real_s(      -9.056250029794158606932797 ) * liVS[42] +
                   li_Real_s(      10.062500048821634734963482 ) * liVS[43] +
                   li_Real_s(      -7.546875047463875496589480 ) * liVS[44] +
                   li_Real_s(       3.622500027248760190445864 ) * liVS[45] +
                   li_Real_s(      -1.006250008542771912090075 ) * liVS[46] +
                   li_Real_s(       0.123214286844747000770894 ) * liVS[47] +
                   li_Real_s(       0.684226190947875956283042 ) * liVS[48] +
                   li_Real_s(      -1.847222225457541000537276 ) * liVS[49] +
                   li_Real_s(       2.770833342907319085668405 ) * liVS[50] +
                   li_Real_s(      -3.078703719381994829973337 ) * liVS[51] +
                   li_Real_s(       2.309027793005110851254358 ) * liVS[52] +
                   li_Real_s(      -1.108333342066313065288341 ) * liVS[53] +
                   li_Real_s(       0.307870373105676442548884 ) * liVS[54] +
                   li_Real_s(      -0.037698413060133439955735 ) * liVS[55] +
                   li_Real_s(      -0.090029761967963750635136 ) * liVS[56] +
                   li_Real_s(       0.243055555990299460011883 ) * liVS[57] +
                   li_Real_s(      -0.364583334621044086176767 ) * liVS[58] +
                   li_Real_s(       0.405092594701361718989574 ) * liVS[59] +
                   li_Real_s(      -0.303819446491643052610243 ) * liVS[60] +
                   li_Real_s(       0.145833334506537459773767 ) * liVS[61] +
                   li_Real_s(      -0.040509259626368976370259 ) * liVS[62] +
                   li_Real_s(       0.004960317508818867793252 ) * liVS[63];
        liVC[14] = li_Real_s(      -0.014404761910646435296712 ) * liVS[ 0] +
                   li_Real_s(       0.038888888924385967005293 ) * liVS[ 1] +
                   li_Real_s(      -0.058333333430849981793287 ) * liVS[ 2] +
                   li_Real_s(       0.064814814968263245731350 ) * liVS[ 3] +
                   li_Real_s(      -0.048611111256405825642446 ) * liVS[ 4] +
                   li_Real_s(       0.023333333414600088140389 ) * liVS[ 5] +
                   li_Real_s(      -0.006481481506091291144855 ) * liVS[ 6] +
                   li_Real_s(       0.000793650796744566067176 ) * liVS[ 7] +
                   li_Real_s(       0.097232142915668118288153 ) * liVS[ 8] +
                   li_Real_s(      -0.262500000377247744154374 ) * liVS[ 9] +
                   li_Real_s(       0.393750001075369615577415 ) * liVS[10] +
                   li_Real_s(      -0.437500001722711751028783 ) * liVS[11] +
                   li_Real_s(       0.328125001649731407038502 ) * liVS[12] +
                   li_Real_s(      -0.157500000935730688489400 ) * liVS[13] +
                   li_Real_s(       0.043750000289854984458771 ) * liVS[14] +
                   li_Real_s(      -0.005357142894933941690283 ) * liVS[15] +
                   li_Real_s(      -0.280892857349566482660919 ) * liVS[16] +
                   li_Real_s(       0.758333334695539873493431 ) * liVS[17] +
                   li_Real_s(      -1.137500003924288183299041 ) * liVS[18] +
                   li_Real_s(       1.263888895202114337479316 ) * liVS[19] +
                   li_Real_s(      -0.947916672726471842835849 ) * liVS[20] +
                   li_Real_s(       0.455000003448395862193365 ) * liVS[21] +
                   li_Real_s(      -0.126388889963381600978209 ) * liVS[22] +
                   li_Real_s(       0.015476190617657925585604 ) * liVS[23] +
                   li_Real_s(       0.450148809890471790495781 ) * liVS[24] +
                   li_Real_s(      -1.215277780217707270793426 ) * liVS[25] +
                   li_Real_s(       1.822916673724519487365114 ) * liVS[26] +
                   li_Real_s(      -2.025462974331264209837400 ) * liVS[27] +
                   li_Real_s(       1.519097233138579694511350 ) * liVS[28] +
                   li_Real_s(      -0.729166672883657551729186 ) * liVS[29] +
                   li_Real_s(       0.202546298237046895618363 ) * liVS[30] +
                   li_Real_s(      -0.024801587557989366455979 ) * liVS[31] +
                   li_Real_s(      -0.432142857511714550966531 ) * liVS[32] +
                   li_Real_s(       1.166666669133562983518004 ) * liVS[33] +
                   li_Real_s(      -1.750000007148266201184583 ) * liVS[34] +
                   li_Real_s(       1.944444455960068429689613 ) * liVS[35] +
                   li_Real_s(      -1.458333344387532370944882 ) * liVS[36] +
                   li_Real_s(       0.700000006294234644599328 ) * liVS[37] +
                   li_Real_s(      -0.194444446409782034645986 ) * liVS[38] +
                   li_Real_s(       0.023809524069428711356977 ) * liVS[39] +
                   li_Real_s(       0.248482143072332184452478 ) * liVS[40] +
                   li_Real_s(      -0.670833334776932077581080 ) * liVS[41] +
                   li_Real_s(       1.006250004186703073116860 ) * liVS[42] +
                   li_Real_s(      -1.118055562298933436338189 ) * liVS[43] +
                   li_Real_s(       0.838541673136241283437187 ) * liVS[44] +
                   li_Real_s(      -0.402500003681566076441101 ) * liVS[45] +
                   li_Real_s(       0.111805556704558328728893 ) * liVS[46] +
                   li_Real_s(      -0.013690476342402837020562 ) * liVS[47] +
                   li_Real_s(      -0.079226190544344377020280 ) * liVS[48] +
                   li_Real_s(       0.213888889347179700362034 ) * liVS[49] +
                   li_Real_s(      -0.320833334663410130183081 ) * liVS[50] +
                   li_Real_s(       0.356481483623580341912884 ) * liVS[51] +
                   li_Real_s(      -0.267361113165280439574190 ) * liVS[52] +
                   li_Real_s(       0.128333334501523232962228 ) * liVS[53] +
                   li_Real_s(      -0.035648148512468669224518 ) * liVS[54] +
                   li_Real_s(       0.004365079413220174231469 ) * liVS[55] +
                   li_Real_s(       0.010803571437668191279613 ) * liVS[56] +
                   li_Real_s(      -0.029166666727958721894165 ) * liVS[57] +
                   li_Real_s(       0.043750000178016723584307 ) * liVS[58] +
                   li_Real_s(      -0.048611111397832695857346 ) * liVS[59] +
                   li_Real_s(       0.036458333608205328868479 ) * liVS[60] +
                   li_Real_s(      -0.017500000156226165615170 ) * liVS[61] +
                   li_Real_s(       0.004861111159792555480585 ) * liVS[62] +
                   li_Real_s(      -0.000595238101664952168335 ) * liVS[63];
        liVC[15] = li_Real_s(       0.000514455782690853569150 ) * liVS[ 0] +
                   li_Real_s(      -0.001388888891204143050118 ) * liVS[ 1] +
                   li_Real_s(       0.002083333339664362321431 ) * liVS[ 2] +
                   li_Real_s(      -0.002314814824640971935299 ) * liVS[ 3] +
                   li_Real_s(       0.001736111120295150828108 ) * liVS[ 4] +
                   li_Real_s(      -0.000833333338438491111333 ) * liVS[ 5] +
                   li_Real_s(       0.000231481483033250511755 ) * liVS[ 6] +
                   li_Real_s(      -0.000028344671400029781971 ) * liVS[ 7] +
                   li_Real_s(      -0.003601190479460048693428 ) * liVS[ 8] +
                   li_Real_s(       0.009722222243183193182703 ) * liVS[ 9] +
                   li_Real_s(      -0.014583333392384526319696 ) * liVS[10] +
                   li_Real_s(       0.016203703797086368987301 ) * liVS[11] +
                   li_Real_s(      -0.012152777866255593322564 ) * liVS[12] +
                   li_Real_s(       0.005833333383182293628566 ) * liVS[13] +
                   li_Real_s(      -0.001620370385778208306704 ) * liVS[14] +
                   li_Real_s(       0.000198412700426531252162 ) * liVS[15] +
                   li_Real_s(       0.010803571439440703971790 ) * liVS[16] +
                   li_Real_s(      -0.029166666737661453567831 ) * liVS[17] +
                   li_Real_s(       0.043750000202218738976079 ) * liVS[18] +
                   li_Real_s(      -0.048611111432918019392702 ) * liVS[19] +
                   li_Real_s(       0.036458333639560927541545 ) * liVS[20] +
                   li_Real_s(      -0.017500000173272176051675 ) * liVS[21] +
                   li_Real_s(       0.004861111164961656738726 ) * liVS[22] +
                   li_Real_s(      -0.000595238102330399032613 ) * liVS[23] +
                   li_Real_s(      -0.018005952399660470064635 ) * liVS[24] +
                   li_Real_s(       0.048611111234394155378169 ) * liVS[25] +
                   li_Real_s(      -0.072916667019487818590662 ) * liVS[26] +
                   li_Real_s(       0.081018519081342207921459 ) * liVS[27] +
                   li_Real_s(      -0.060763889425251405018713 ) * liVS[28] +
                   li_Real_s(       0.029166666970584981255499 ) * liVS[29] +
                   li_Real_s(      -0.008101851946484196764775 ) * liVS[30] +
                   li_Real_s(       0.000992063504562550220467 ) * liVS[31] +
                   li_Real_s(       0.018005952399463093227538 ) * liVS[32] +
                   li_Real_s(      -0.048611111233657189334423 ) * liVS[33] +
                   li_Real_s(       0.072916667018117442555791 ) * liVS[34] +
                   li_Real_s(      -0.081018519079590317621964 ) * liVS[35] +
                   li_Real_s(       0.060763889423709832593445 ) * liVS[36] +
                   li_Real_s(      -0.029166666969757240601702 ) * liVS[37] +
                   li_Real_s(       0.008101851946256172568028 ) * liVS[38] +
                   li_Real_s(      -0.000992063504541819407567 ) * liVS[39] +
                   li_Real_s(      -0.010803571439265479553482 ) * liVS[40] +
                   li_Real_s(       0.029166666737656520014266 ) * liVS[41] +
                   li_Real_s(      -0.043750000203796504671949 ) * liVS[42] +
                   li_Real_s(       0.048611111436502492388900 ) * liVS[43] +
                   li_Real_s(      -0.036458333643421762804149 ) * liVS[44] +
                   li_Real_s(       0.017500000175661414164585 ) * liVS[45] +
                   li_Real_s(      -0.004861111165801975872447 ) * liVS[46] +
                   li_Real_s(       0.000595238102465305007893 ) * liVS[47] +
                   li_Real_s(       0.003601190479557777809894 ) * liVS[48] +
                   li_Real_s(      -0.009722222244618841657804 ) * liVS[49] +
                   li_Real_s(       0.014583333397672638281906 ) * liVS[50] +
                   li_Real_s(      -0.016203703806432823153472 ) * liVS[51] +
                   li_Real_s(       0.012152777875643901162039 ) * liVS[52] +
                   li_Real_s(      -0.005833333388743242031715 ) * liVS[53] +
                   li_Real_s(       0.001620370387609808282559 ) * liVS[54] +
                   li_Real_s(      -0.000198412700689227367024 ) * liVS[55] +
                   li_Real_s(      -0.000514455782761054358776 ) * liVS[56] +
                   li_Real_s(       0.001388888891873716821546 ) * liVS[57] +
                   li_Real_s(      -0.002083333341912776449922 ) * liVS[58] +
                   li_Real_s(       0.002314814828514245165225 ) * liVS[59] +
                   li_Real_s(      -0.001736111124158534399498 ) * liVS[60] +
                   li_Real_s(       0.000833333340716454194254 ) * liVS[61] +
                   li_Real_s(      -0.000231481483776690977194 ) * liVS[62] +
                   li_Real_s(       0.000028344671504622331870 ) * liVS[63];
        liVC[16] = li_Real_s(       2.605555556137634010127613 ) * liVS[ 0] +
                   li_Real_s(     -11.150000003808660409276854 ) * liVS[ 1] +
                   li_Real_s(      21.975000010768475533495803 ) * liVS[ 2] +
                   li_Real_s(     -26.361111128064521835767664 ) * liVS[ 3] +
                   li_Real_s(      20.500000016006751479835657 ) * liVS[ 4] +
                   li_Real_s(     -10.050000009058793537519705 ) * liVS[ 5] +
                   li_Real_s(       2.830555558402628513192667 ) * liVS[ 6] +
                   li_Real_s(      -0.350000000383515030843995 ) * liVS[ 7] +
                   li_Real_s(      -0.000000003670194978780301 ) * liVS[ 8] +
                   li_Real_s(       0.000000024088454276037542 ) * liVS[ 9] +
                   li_Real_s(      -0.000000068201969487242849 ) * liVS[10] +
                   li_Real_s(       0.000000107402172016630605 ) * liVS[11] +
                   li_Real_s(      -0.000000101361670889663439 ) * liVS[12] +
                   li_Real_s(       0.000000057319958311661445 ) * liVS[13] +
                   li_Real_s(      -0.000000017999406615354459 ) * liVS[14] +
                   li_Real_s(       0.000000002422657366711421 ) * liVS[15] +
                   li_Real_s(       0.000000009884682725752305 ) * liVS[16] +
                   li_Real_s(      -0.000000065086909435010333 ) * liVS[17] +
                   li_Real_s(       0.000000184573783031190906 ) * liVS[18] +
                   li_Real_s(      -0.000000290778075481638305 ) * liVS[19] +
                   li_Real_s(       0.000000274329271133572276 ) * liVS[20] +
                   li_Real_s(      -0.000000155017625243215858 ) * liVS[21] +
                   li_Real_s(       0.000000048635344404723104 ) * liVS[22] +
                   li_Real_s(      -0.000000006540471135374068 ) * liVS[23] +
                   li_Real_s(      -0.000000014696790781823158 ) * liVS[24] +
                   li_Real_s(       0.000000097118268872987125 ) * liVS[25] +
                   li_Real_s(      -0.000000275941996606300052 ) * liVS[26] +
                   li_Real_s(       0.000000435031612787514823 ) * liVS[27] +
                   li_Real_s(      -0.000000410382216569854218 ) * liVS[28] +
                   li_Real_s(       0.000000231767721298294413 ) * liVS[29] +
                   li_Real_s(      -0.000000072660149876610571 ) * liVS[30] +
                   li_Real_s(       0.000000009763550875791605 ) * liVS[31] +
                   li_Real_s(       0.000000013032967185502724 ) * liVS[32] +
                   li_Real_s(      -0.000000086460865174670233 ) * liVS[33] +
                   li_Real_s(       0.000000246233630955157697 ) * liVS[34] +
                   li_Real_s(      -0.000000388616228942759129 ) * liVS[35] +
                   li_Real_s(       0.000000366674068489059592 ) * liVS[36] +
                   li_Real_s(      -0.000000207017489448666477 ) * liVS[37] +
                   li_Real_s(       0.000000064864132945010512 ) * liVS[38] +
                   li_Real_s(      -0.000000008710216008634699 ) * liVS[39] +
                   li_Real_s(      -0.000000006908319296576560 ) * liVS[40] +
                   li_Real_s(       0.000000046022092969421475 ) * liVS[41] +
                   li_Real_s(      -0.000000131416631466340371 ) * liVS[42] +
                   li_Real_s(       0.000000207699722041154986 ) * liVS[43] +
                   li_Real_s(      -0.000000196069788458400553 ) * liVS[44] +
                   li_Real_s(       0.000000110688427823317976 ) * liVS[45] +
                   li_Real_s(      -0.000000034668515098230032 ) * liVS[46] +
                   li_Real_s(       0.000000004653011485653155 ) * liVS[47] +
                   li_Real_s(       0.000000002033177636907905 ) * liVS[48] +
                   li_Real_s(      -0.000000013602840738909941 ) * liVS[49] +
                   li_Real_s(       0.000000038953638385455872 ) * liVS[50] +
                   li_Real_s(      -0.000000061664340837799843 ) * liVS[51] +
                   li_Real_s(       0.000000058251729151927490 ) * liVS[52] +
                   li_Real_s(      -0.000000032888111111069646 ) * liVS[53] +
                   li_Real_s(       0.000000010298371411492682 ) * liVS[54] +
                   li_Real_s(      -0.000000001381623898004546 ) * liVS[55] +
                   li_Real_s(      -0.000000000257310143164569 ) * liVS[56] +
                   li_Real_s(       0.000000001728628596132874 ) * liVS[57] +
                   li_Real_s(      -0.000000004964099702512304 ) * liVS[58] +
                   li_Real_s(       0.000000007871238504774880 ) * liVS[59] +
                   li_Real_s(      -0.000000007441276225682825 ) * liVS[60] +
                   li_Real_s(       0.000000004201942257521208 ) * liVS[61] +
                   li_Real_s(      -0.000000001315551835469501 ) * liVS[62] +
                   li_Real_s(       0.000000000176428548400264 ) * liVS[63];
        liVC[17] = li_Real_s(      -6.755833333644943650142523 ) * liVS[ 0] +
                   li_Real_s(      28.910357144812955709767266 ) * liVS[ 1] +
                   li_Real_s(     -56.978035720471972069844924 ) * liVS[ 2] +
                   li_Real_s(      68.350595249412947396194795 ) * liVS[ 3] +
                   li_Real_s(     -53.153571440714898699297919 ) * liVS[ 4] +
                   li_Real_s(      26.058214293210852474658168 ) * liVS[ 5] +
                   li_Real_s(      -7.339226192937208459454723 ) * liVS[ 6] +
                   li_Real_s(       0.907500000331729950175941 ) * liVS[ 7] +
                   li_Real_s(      18.238888893970994331539259 ) * liVS[ 8] +
                   li_Real_s(     -78.050000034173720564467658 ) * liVS[ 9] +
                   li_Real_s(     153.825000103184635236175382 ) * liVS[10] +
                   li_Real_s(    -184.527777952375373615723220 ) * liVS[11] +
                   li_Real_s(     143.500000175142133684857981 ) * liVS[12] +
                   li_Real_s(     -70.350000103300317277899012 ) * liVS[13] +
                   li_Real_s(      19.813888922006775317186111 ) * liVS[14] +
                   li_Real_s(      -2.450000004455148427950917 ) * liVS[15] +
                   li_Real_s(     -27.358333353956467703937960 ) * liVS[16] +
                   li_Real_s(     117.075000140288096872609458 ) * liVS[17] +
                   li_Real_s(    -230.737500419027497855495312 ) * liVS[18] +
                   li_Real_s(     276.791667363121177913853899 ) * liVS[19] +
                   li_Real_s(    -215.250000687065210058790399 ) * liVS[20] +
                   li_Real_s(     105.525000400269050260249060 ) * liVS[21] +
                   li_Real_s(     -29.720833460742724696501682 ) * liVS[22] +
                   li_Real_s(       3.675000017113603689722368 ) * liVS[23] +
                   li_Real_s(      30.398148186823817695767502 ) * liVS[24] +
                   li_Real_s(    -130.083333597835689943167381 ) * liVS[25] +
                   li_Real_s(     256.375000786920338669006014 ) * liVS[26] +
                   li_Real_s(    -307.546297594884094905864913 ) * liVS[27] +
                   li_Real_s(     239.166667938784570424104459 ) * liVS[28] +
                   li_Real_s(    -117.250000737039002274286759 ) * liVS[29] +
                   li_Real_s(      33.023148381925182093254989 ) * liVS[30] +
                   li_Real_s(      -4.083333364694766487446032 ) * liVS[31] +
                   li_Real_s(     -22.798611151044099187856773 ) * liVS[32] +
                   li_Real_s(      97.562500274138415079505648 ) * liVS[33] +
                   li_Real_s(    -192.281250814331343690355425 ) * liVS[34] +
                   li_Real_s(     230.659723561287876236747252 ) * liVS[35] +
                   li_Real_s(    -179.375001306864959360609646 ) * liVS[36] +
                   li_Real_s(      87.937500754767100374920119 ) * liVS[37] +
                   li_Real_s(     -24.767361349953795013334457 ) * liVS[38] +
                   li_Real_s(       3.062500032000812666410638 ) * liVS[39] +
                   li_Real_s(      10.943333356960991409323469 ) * liVS[40] +
                   li_Real_s(     -46.830000162726904022747476 ) * liVS[41] +
                   li_Real_s(      92.295000483178540662265732 ) * liVS[42] +
                   li_Real_s(    -110.716667459661991301800299 ) * liVS[43] +
                   li_Real_s(      86.100000772158495010444312 ) * liVS[44] +
                   li_Real_s(     -42.210000444997021418203076 ) * liVS[45] +
                   li_Real_s(      11.888333473897269243479968 ) * liVS[46] +
                   li_Real_s(      -1.470000018809335617930856 ) * liVS[47] +
                   li_Real_s(      -3.039814822370757951830456 ) * liVS[48] +
                   li_Real_s(      13.008333385529212478104455 ) * liVS[49] +
                   li_Real_s(     -25.637500155020006786799058 ) * liVS[50] +
                   li_Real_s(      30.754629883791807287707343 ) * liVS[51] +
                   li_Real_s(     -23.916666913790209036960732 ) * liVS[52] +
                   li_Real_s(      11.725000142197885111272626 ) * liVS[53] +
                   li_Real_s(      -3.302314859665486324047379 ) * liVS[54] +
                   li_Real_s(       0.408333339327541011698486 ) * liVS[55] +
                   li_Real_s(       0.372222223238424021474202 ) * liVS[56] +
                   li_Real_s(      -1.592857149896806934208371 ) * liVS[57] +
                   li_Real_s(       3.139285735207153038572869 ) * liVS[58] +
                   li_Real_s(      -3.765873050160720936219150 ) * liVS[59] +
                   li_Real_s(       2.928571461879393211802380 ) * liVS[60] +
                   li_Real_s(      -1.435714304858241252560447 ) * liVS[61] +
                   li_Real_s(       0.404365085395749446206537 ) * liVS[62] +
                   li_Real_s(      -0.050000000804942601462244 ) * liVS[63];
        liVC[18] = li_Real_s(       6.788919751358434950816445 ) * liVS[ 0] +
                   li_Real_s(     -29.051944433223013675160473 ) * liVS[ 1] +
                   li_Real_s(      57.257083303858877343373024 ) * liVS[ 2] +
                   li_Real_s(     -68.685339464440801293676486 ) * liVS[ 3] +
                   li_Real_s(      53.413888853643200604892627 ) * liVS[ 4] +
                   li_Real_s(     -26.185833315164586565515492 ) * liVS[ 5] +
                   li_Real_s(       7.375169747691979083015212 ) * liVS[ 6] +
                   li_Real_s(      -0.911944443724571840448334 ) * liVS[ 7] +
                   li_Real_s(     -29.051944440401882729929639 ) * liVS[ 8] +
                   li_Real_s(     124.322499977148183347708255 ) * liVS[ 9] +
                   li_Real_s(    -245.021249955009011500806082 ) * liVS[10] +
                   li_Real_s(     293.926388852689171926613199 ) * liVS[11] +
                   li_Real_s(    -228.574999994876264963750145 ) * liVS[12] +
                   li_Real_s(     112.057500009093445214602980 ) * liVS[13] +
                   li_Real_s(     -31.560694449336281763862644 ) * liVS[14] +
                   li_Real_s(       3.902500000692668891133508 ) * liVS[15] +
                   li_Real_s(      57.257083337400786149373744 ) * liVS[16] +
                   li_Real_s(    -245.021250042376976807645406 ) * liVS[17] +
                   li_Real_s(     482.900625183033639586938079 ) * liVS[18] +
                   li_Real_s(    -579.285417057750009917072020 ) * liVS[19] +
                   li_Real_s(     450.487500452460210453864420 ) * liVS[20] +
                   li_Real_s(    -220.848750289364147647575010 ) * liVS[21] +
                   li_Real_s(      62.201458429597082044892886 ) * liVS[22] +
                   li_Real_s(      -7.691250013000683338759700 ) * liVS[23] +
                   li_Real_s(     -68.685339529013390347245149 ) * liVS[24] +
                   li_Real_s(     293.926389068905962176359026 ) * liVS[25] +
                   li_Real_s(    -579.285417281996274141420145 ) * liVS[26] +
                   li_Real_s(     694.908180145823052953346632 ) * liVS[27] +
                   li_Real_s(    -540.402778977642356039723381 ) * liVS[28] +
                   li_Real_s(     264.929167396167713377508335 ) * liVS[29] +
                   li_Real_s(     -74.616589742974042565037962 ) * liVS[30] +
                   li_Real_s(       9.226388920729258202868550 ) * liVS[31] +
                   li_Real_s(      53.413888920485817379812943 ) * liVS[32] +
                   li_Real_s(    -228.575000239186380213141092 ) * liVS[33] +
                   li_Real_s(     450.487500779990114097017795 ) * liVS[34] +
                   li_Real_s(    -540.402779160878822040103842 ) * liVS[35] +
                   li_Real_s(     420.250001425343668870482361 ) * liVS[36] +
                   li_Real_s(    -206.025000851842690963167115 ) * liVS[37] +
                   li_Real_s(      58.026389162818730937942746 ) * liVS[38] +
                   li_Real_s(      -7.175000036730523333972087 ) * liVS[39] +
                   li_Real_s(     -26.185833354783198956283741 ) * liVS[40] +
                   li_Real_s(     112.057500160083193918580946 ) * liVS[41] +
                   li_Real_s(    -220.848750512147120161898783 ) * liVS[42] +
                   li_Real_s(     264.929167559614711535687093 ) * liVS[43] +
                   li_Real_s(    -206.025000908535247390318546 ) * liVS[44] +
                   li_Real_s(     101.002500538211506864172406 ) * liVS[45] +
                   li_Real_s(     -28.447083505469294095746591 ) * liVS[46] +
                   li_Real_s(       3.517500023025572186696763 ) * liVS[47] +
                   li_Real_s(       7.375169760476069313881453 ) * liVS[48] +
                   li_Real_s(     -31.560694499297753168320924 ) * liVS[49] +
                   li_Real_s(      62.201458507219001603516517 ) * liVS[50] +
                   li_Real_s(     -74.616589806687983355004690 ) * liVS[51] +
                   li_Real_s(      58.026389192479655321221799 ) * liVS[52] +
                   li_Real_s(     -28.447083512226001289491251 ) * liVS[53] +
                   li_Real_s(       8.012044810089264501584694 ) * liVS[54] +
                   li_Real_s(      -0.990694452052281349097029 ) * liVS[55] +
                   li_Real_s(      -0.911944445476990495080827 ) * liVS[56] +
                   li_Real_s(       3.902500007659380543145744 ) * liVS[57] +
                   li_Real_s(      -7.691250024182892275348422 ) * liVS[58] +
                   li_Real_s(       9.226388930495886597782373 ) * liVS[59] +
                   li_Real_s(      -7.175000041865075672831153 ) * liVS[60] +
                   li_Real_s(       3.517500024587151497001969 ) * liVS[61] +
                   li_Real_s(      -0.990694452257528723748692 ) * liVS[62] +
                   li_Real_s(       0.122500001040236838889541 ) * liVS[63];
        liVC[19] = li_Real_s(      -3.499405862758408147783484 ) * liVS[ 0] +
                   li_Real_s(      14.975069435185048405401176 ) * liVS[ 1] +
                   li_Real_s(     -29.513645809150517607122310 ) * liVS[ 2] +
                   li_Real_s(      35.404436694284015629818896 ) * liVS[ 3] +
                   li_Real_s(     -27.532638860131456226554292 ) * liVS[ 4] +
                   li_Real_s(      13.497708318491618229018059 ) * liVS[ 5] +
                   li_Real_s(      -3.801593359768846269730602 ) * liVS[ 6] +
                   li_Real_s(       0.470069443848657897433441 ) * liVS[ 7] +
                   li_Real_s(      18.470493823623854723336990 ) * liVS[ 8] +
                   li_Real_s(     -79.041111091611583105986938 ) * liVS[ 9] +
                   li_Real_s(     155.778333295588254259200767 ) * liVS[10] +
                   li_Real_s(    -186.870987624655469971912680 ) * liVS[11] +
                   li_Real_s(     145.322222218759634415619075 ) * liVS[12] +
                   li_Real_s(     -71.243333341280035142517590 ) * liVS[13] +
                   li_Real_s(      20.065493831239841426850035 ) * liVS[14] +
                   li_Real_s(      -2.481111111664517920871731 ) * liVS[15] +
                   li_Real_s(     -42.655115743251201365637826 ) * liVS[16] +
                   li_Real_s(     182.534791698078038280073088 ) * liVS[17] +
                   li_Real_s(    -359.749062644774312502704561 ) * liVS[18] +
                   li_Real_s(     431.553356798471327238075901 ) * liVS[19] +
                   li_Real_s(    -335.602083703691448590689106 ) * liVS[20] +
                   li_Real_s(     164.526875237528855677737738 ) * liVS[21] +
                   li_Real_s(     -46.338553319628360327442351 ) * liVS[22] +
                   li_Real_s(       5.729791677267087379732402 ) * liVS[23] +
                   li_Real_s(      56.308950634476218510826584 ) * liVS[24] +
                   li_Real_s(    -240.963889029615131676109741 ) * liVS[25] +
                   li_Real_s(     474.904167158348059274430852 ) * liVS[26] +
                   li_Real_s(    -569.692902151492489792872220 ) * liVS[27] +
                   li_Real_s(     443.027778754807002314919373 ) * liVS[28] +
                   li_Real_s(    -217.191667262248330416696263 ) * liVS[29] +
                   li_Real_s(      61.171450810520980212459108 ) * liVS[30] +
                   li_Real_s(      -7.563888914796235596327278 ) * liVS[31] +
                   li_Real_s(     -46.049575641567969341849675 ) * liVS[32] +
                   li_Real_s(     197.060764077186689746667980 ) * liVS[33] +
                   li_Real_s(    -388.377604790054135719401529 ) * liVS[34] +
                   li_Real_s(     465.896027350384883902734146 ) * liVS[35] +
                   li_Real_s(    -362.309028933957733897841536 ) * liVS[36] +
                   li_Real_s(     177.619792359284701888100244 ) * liVS[37] +
                   li_Real_s(     -50.026138339951351952095138 ) * liVS[38] +
                   li_Real_s(       6.185763918674808792275144 ) * liVS[39] +
                   li_Real_s(      23.276296312876013416826027 ) * liVS[40] +
                   li_Real_s(     -99.606666792790207409780123 ) * liVS[41] +
                   li_Real_s(     196.310000408549171879712958 ) * liVS[42] +
                   li_Real_s(    -235.492593310690750740832300 ) * liVS[43] +
                   li_Real_s(     183.133334067536964084865758 ) * liVS[44] +
                   li_Real_s(     -89.780000435876857522998762 ) * liVS[45] +
                   li_Real_s(      25.286296435661341774903121 ) * liVS[46] +
                   li_Real_s(      -3.126666685265561795858957 ) * liVS[47] +
                   li_Real_s(      -6.691211425464459239265125 ) * liVS[48] +
                   li_Real_s(      28.633819487569770245727341 ) * liVS[49] +
                   li_Real_s(     -56.433020971591560055458103 ) * liVS[50] +
                   li_Real_s(      67.696798080250005114066880 ) * liVS[51] +
                   li_Real_s(     -52.645139133221846350352280 ) * liVS[52] +
                   li_Real_s(      25.808958477600228320625320 ) * liVS[53] +
                   li_Real_s(      -7.269023965704974443724495 ) * liVS[54] +
                   li_Real_s(       0.898819450562768906820565 ) * liVS[55] +
                   li_Real_s(       0.839567902030125878809486 ) * liVS[56] +
                   li_Real_s(      -3.592777783774987909737320 ) * liVS[57] +
                   li_Real_s(       7.080833352476957998078433 ) * liVS[58] +
                   li_Real_s(      -8.494135835650098442783928 ) * liVS[59] +
                   li_Real_s(       6.605555589095942536914663 ) * liVS[60] +
                   li_Real_s(      -3.238333353070526499095649 ) * liVS[61] +
                   li_Real_s(       0.912067907503384844858374 ) * liVS[62] +
                   li_Real_s(      -0.112777778610677614778979 ) * liVS[63];
        liVC[20] = li_Real_s(       1.013271604648622314925888 ) * liVS[ 0] +
                   li_Real_s(      -4.336111109270461838605115 ) * liVS[ 1] +
                   li_Real_s(       8.545833328845382936833630 ) * liVS[ 2] +
                   li_Real_s(     -10.251543204237890449803672 ) * liVS[ 3] +
                   li_Real_s(       7.972222218126951531758095 ) * liVS[ 4] +
                   li_Real_s(      -3.908333331504053376193042 ) * liVS[ 5] +
                   li_Real_s(       1.100771604430625494330798 ) * liVS[ 6] +
                   li_Real_s(      -0.136111111039316057258475 ) * liVS[ 7] +
                   li_Real_s(      -6.025347222559283011378284 ) * liVS[ 8] +
                   li_Real_s(      25.784375003666021086701221 ) * liVS[ 9] +
                   li_Real_s(     -50.817187516955556247921777 ) * liVS[10] +
                   li_Real_s(      60.960069482266419527149992 ) * liVS[11] +
                   li_Real_s(     -47.406250044751679695309576 ) * liVS[12] +
                   li_Real_s(      23.240625028834802634492007 ) * liVS[13] +
                   li_Real_s(      -6.545659731766594546797933 ) * liVS[14] +
                   li_Real_s(       0.809375001265877358491707 ) * liVS[15] +
                   li_Real_s(      15.416203709364822316274513 ) * liVS[16] +
                   li_Real_s(     -65.970833376471915698857629 ) * liVS[17] +
                   li_Real_s(     130.018750143629375770615297 ) * liVS[18] +
                   li_Real_s(    -155.969907667434085851709824 ) * liVS[19] +
                   li_Real_s(     121.291666938804667097429046 ) * liVS[20] +
                   li_Real_s(     -59.462500164093214038985025 ) * liVS[21] +
                   li_Real_s(      16.747453756593536411401146 ) * liVS[22] +
                   li_Real_s(      -2.070833340393178900740168 ) * liVS[23] +
                   li_Real_s(     -22.056751557042275635467377 ) * liVS[24] +
                   li_Real_s(      94.387847324233533186088607 ) * liVS[25] +
                   li_Real_s(    -186.024479491674867404071847 ) * liVS[26] +
                   li_Real_s(     223.154128652882661754119908 ) * liVS[27] +
                   li_Real_s(    -173.538195021291187458700733 ) * liVS[28] +
                   li_Real_s(      85.076042008619552348136494 ) * liVS[29] +
                   li_Real_s(     -23.961439152547782782676222 ) * liVS[30] +
                   li_Real_s(       2.962847236820358443054602 ) * liVS[31] +
                   li_Real_s(      19.107407423401525647932431 ) * liVS[32] +
                   li_Real_s(     -81.766666783260546935707680 ) * liVS[33] +
                   li_Real_s(     161.150000365198422969115200 ) * liVS[34] +
                   li_Real_s(    -193.314815441397456652339315 ) * liVS[35] +
                   li_Real_s(     150.333333963732059146423126 ) * liVS[36] +
                   li_Real_s(     -73.700000370656084669462871 ) * liVS[37] +
                   li_Real_s(      20.757407525381196933267347 ) * liVS[38] +
                   li_Real_s(      -2.566666682399088017518807 ) * liVS[39] +
                   li_Real_s(     -10.042245380418322042714863 ) * liVS[40] +
                   li_Real_s(      42.973958406210570615257893 ) * liVS[41] +
                   li_Real_s(     -84.695312726319571083877236 ) * liVS[42] +
                   li_Real_s(     101.600116125797100607996981 ) * liVS[43] +
                   li_Real_s(     -79.010417051425847034806793 ) * liVS[44] +
                   li_Real_s(      38.734375225088648164728511 ) * liVS[45] +
                   li_Real_s(     -10.909432941771147085319171 ) * liVS[46] +
                   li_Real_s(       1.348958342838529222973420 ) * liVS[47] +
                   li_Real_s(       2.967438274923075880451506 ) * liVS[48] +
                   li_Real_s(     -12.698611135140268402210495 ) * liVS[49] +
                   li_Real_s(      25.027083407632261469188961 ) * liVS[50] +
                   li_Real_s(     -30.022376669026726858646725 ) * liVS[51] +
                   li_Real_s(      23.347222347422302846098319 ) * liVS[52] +
                   li_Real_s(     -11.445833406329622761177234 ) * liVS[53] +
                   li_Real_s(       3.223688294698547451844206 ) * liVS[54] +
                   li_Real_s(      -0.398611114179566072834859 ) * liVS[55] +
                   li_Real_s(      -0.379976852304238832402916 ) * liVS[56] +
                   li_Real_s(       1.626041669945414547271412 ) * liVS[57] +
                   li_Real_s(      -3.204687510120837856675280 ) * liVS[58] +
                   li_Real_s(       3.844328720801520660188544 ) * liVS[59] +
                   li_Real_s(      -2.989583350306560305398307 ) * liVS[60] +
                   li_Real_s(       1.465625009873416262351498 ) * liVS[61] +
                   li_Real_s(      -0.412789354968651878152741 ) * liVS[62] +
                   li_Real_s(       0.051041667079957164787629 ) * liVS[63];
        liVC[21] = li_Real_s(      -0.166466049411368999244587 ) * liVS[ 0] +
                   li_Real_s(       0.712361111301138549833922 ) * liVS[ 1] +
                   li_Real_s(      -1.403958333985609208127698 ) * liVS[ 2] +
                   li_Real_s(       1.684182100034188778181488 ) * liVS[ 3] +
                   li_Real_s(      -1.309722223628087434121881 ) * liVS[ 4] +
                   li_Real_s(       0.642083334206276390432322 ) * liVS[ 5] +
                   li_Real_s(      -0.180841049663987618600913 ) * liVS[ 6] +
                   li_Real_s(       0.022361111147374934660093 ) * liVS[ 7] +
                   li_Real_s(       1.067554012971511312457551 ) * liVS[ 8] +
                   li_Real_s(      -4.568402782137956030794612 ) * liVS[ 9] +
                   li_Real_s(       9.003645846850734102417846 ) * liVS[10] +
                   li_Real_s(     -10.800733047985698220827544 ) * liVS[11] +
                   li_Real_s(       8.399305579165822877030223 ) * liVS[12] +
                   li_Real_s(      -4.117708347288421499854394 ) * liVS[13] +
                   li_Real_s(       1.159741516791309656753128 ) * liVS[14] +
                   li_Real_s(      -0.143402778367299532646939 ) * liVS[15] +
                   li_Real_s(      -2.931250002682840971601763 ) * liVS[16] +
                   li_Real_s(      12.543750018737858553663500 ) * liVS[17] +
                   li_Real_s(     -24.721875057043853018967638 ) * liVS[18] +
                   li_Real_s(      29.656250096112870551223750 ) * liVS[19] +
                   li_Real_s(     -23.062500095626063512099790 ) * liVS[20] +
                   li_Real_s(      11.306250055869806914188302 ) * liVS[21] +
                   li_Real_s(      -3.184375017728883427992059 ) * liVS[22] +
                   li_Real_s(       0.393750002361102247050439 ) * liVS[23] +
                   li_Real_s(       4.469251548370792193054513 ) * liVS[24] +
                   li_Real_s(     -19.125347258263879268724850 ) * liVS[25] +
                   li_Real_s(      37.693229275536729971918248 ) * liVS[26] +
                   li_Real_s(     -45.216628268121738187801384 ) * liVS[27] +
                   li_Real_s(      35.163194623770046121080668 ) * liVS[28] +
                   li_Real_s(     -17.238541770864241442495768 ) * liVS[29] +
                   li_Real_s(       4.855189076189279973050361 ) * liVS[30] +
                   li_Real_s(      -0.600347226616969043000438 ) * liVS[31] +
                   li_Real_s(      -4.089274696774609196836536 ) * liVS[32] +
                   li_Real_s(      17.499305593366347721939746 ) * liVS[33] +
                   li_Real_s(     -34.488541780367839351129078 ) * liVS[34] +
                   li_Real_s(      41.372299571452401778515195 ) * liVS[35] +
                   li_Real_s(     -32.173611296485972843584022 ) * liVS[36] +
                   li_Real_s(      15.772916773991711991698139 ) * liVS[37] +
                   li_Real_s(      -4.442399725250601782988724 ) * liVS[38] +
                   li_Real_s(       0.549305560068553688779502 ) * liVS[39] +
                   li_Real_s(       2.247291669907180100551614 ) * liVS[40] +
                   li_Real_s(      -9.616875022617506374444929 ) * liVS[41] +
                   li_Real_s(      18.953437567821687537161779 ) * liVS[42] +
                   li_Real_s(     -22.736458445516880999548448 ) * liVS[43] +
                   li_Real_s(      17.681250109824677707592855 ) * liVS[44] +
                   li_Real_s(      -8.668125063411558528514433 ) * liVS[45] +
                   li_Real_s(       2.441354186649705049205750 ) * liVS[46] +
                   li_Real_s(      -0.301875002657293445285092 ) * liVS[47] +
                   li_Real_s(      -0.687577161536301773026025 ) * liVS[48] +
                   li_Real_s(       2.942361118391225716095505 ) * liVS[49] +
                   li_Real_s(      -5.798958355130384489939388 ) * liVS[50] +
                   li_Real_s(       6.956404356962881507797647 ) * liVS[51] +
                   li_Real_s(      -5.409722257362545860814862 ) * liVS[52] +
                   li_Real_s(       2.652083353581026869960624 ) * liVS[53] +
                   li_Real_s(      -0.746952166862612898512452 ) * liVS[54] +
                   li_Real_s(       0.092361111956707375725273 ) * liVS[55] +
                   li_Real_s(       0.090470679152744537532271 ) * liVS[56] +
                   li_Real_s(      -0.387152778759437654620967 ) * liVS[57] +
                   li_Real_s(       0.763020836270856595007217 ) * liVS[58] +
                   li_Real_s(      -0.915316362867100608013970 ) * liVS[59] +
                   li_Real_s(       0.711805560278798044038240 ) * liVS[60] +
                   li_Real_s(      -0.348958336050605610889619 ) * liVS[61] +
                   li_Real_s(       0.098283179865630287963540 ) * liVS[62] +
                   li_Real_s(      -0.012152777890879207234320 ) * liVS[63];
        liVC[22] = li_Real_s(       0.014475308655563345894279 ) * liVS[ 0] +
                   li_Real_s(      -0.061944444532371023548478 ) * liVS[ 1] +
                   li_Real_s(       0.122083333589210529801505 ) * liVS[ 2] +
                   li_Real_s(      -0.146450617704836805188506 ) * liVS[ 3] +
                   li_Real_s(       0.113888889302051543594985 ) * liVS[ 4] +
                   li_Real_s(      -0.055833333572365528851833 ) * liVS[ 5] +
                   li_Real_s(       0.015725308717041461270547 ) * liVS[ 6] +
                   li_Real_s(      -0.001944444454297145075117 ) * liVS[ 7] +
                   li_Real_s(      -0.097708333465742591705805 ) * liVS[ 8] +
                   li_Real_s(       0.418125000886722775206295 ) * liVS[ 9] +
                   li_Real_s(      -0.824062502604051116961159 ) * liVS[10] +
                   li_Real_s(       0.988541670934666605319308 ) * liVS[11] +
                   li_Real_s(      -0.768750004163832656445265 ) * liVS[12] +
                   li_Real_s(       0.376875002401602654522605 ) * liVS[13] +
                   li_Real_s(      -0.106145834089734714122955 ) * liVS[14] +
                   li_Real_s(       0.013125000100369543787338 ) * liVS[15] +
                   li_Real_s(       0.282268518981726224126305 ) * liVS[16] +
                   li_Real_s(      -1.207916669802405529310363 ) * liVS[17] +
                   li_Real_s(       2.380625009219843768448754 ) * liVS[18] +
                   li_Real_s(      -2.855787052101962153471959 ) * liVS[19] +
                   li_Real_s(       2.220833347976436833448588 ) * liVS[20] +
                   li_Real_s(      -1.088750008425991344651607 ) * liVS[21] +
                   li_Real_s(       0.306643521172183408651790 ) * liVS[22] +
                   li_Real_s(      -0.037916667019831651330719 ) * liVS[23] +
                   li_Real_s(      -0.452353395880046083021853 ) * liVS[24] +
                   li_Real_s(       1.935763894452172229065923 ) * liVS[25] +
                   li_Real_s(      -3.815104183017036465486171 ) * liVS[26] +
                   li_Real_s(       4.576581816776426947512846 ) * liVS[27] +
                   li_Real_s(      -3.559027803617746954500944 ) * liVS[28] +
                   li_Real_s(       1.744791681506844938098766 ) * liVS[29] +
                   li_Real_s(      -0.491415899731253214000049 ) * liVS[30] +
                   li_Real_s(       0.060763889510637172919338 ) * liVS[31] +
                   li_Real_s(       0.434259260081788056595542 ) * liVS[32] +
                   li_Real_s(      -1.858333338935395095958825 ) * liVS[33] +
                   li_Real_s(       3.662500016448610562491695 ) * liVS[34] +
                   li_Real_s(      -4.393518545277391496028940 ) * liVS[35] +
                   li_Real_s(       3.416666692553943729393495 ) * liVS[36] +
                   li_Real_s(      -1.675000014840304363161749 ) * liVS[37] +
                   li_Real_s(       0.471759263922620186715307 ) * liVS[38] +
                   li_Real_s(      -0.058333333953872301691490 ) * liVS[39] +
                   li_Real_s(      -0.249699074554484035104451 ) * liVS[40] +
                   li_Real_s(       1.068541669941780014596588 ) * liVS[41] +
                   li_Real_s(      -2.105937509607119118015817 ) * liVS[42] +
                   li_Real_s(       2.526273163750900607738004 ) * liVS[43] +
                   li_Real_s(      -1.964583348400688489476806 ) * liVS[44] +
                   li_Real_s(       0.963125008623026079845886 ) * liVS[45] +
                   li_Real_s(      -0.271261576779698199146651 ) * liVS[46] +
                   li_Real_s(       0.033541667026283736308123 ) * liVS[47] +
                   li_Real_s(       0.079614197683375409475559 ) * liVS[48] +
                   li_Real_s(      -0.340694445485093422298917 ) * liVS[49] +
                   li_Real_s(       0.671458336384096465820903 ) * liVS[50] +
                   li_Real_s(      -0.805478400010450013013497 ) * liVS[51] +
                   li_Real_s(       0.626388893661162793335961 ) * liVS[52] +
                   li_Real_s(      -0.307083336060568368175439 ) * liVS[53] +
                   li_Real_s(       0.086489198385343746267040 ) * liVS[54] +
                   li_Real_s(      -0.010694444557866722433914 ) * liVS[55] +
                   li_Real_s(      -0.010856481501896109165273 ) * liVS[56] +
                   li_Real_s(       0.046458333472791928098200 ) * liVS[57] +
                   li_Real_s(      -0.091562500408734481815998 ) * liVS[58] +
                   li_Real_s(       0.109837963625471601858408 ) * liVS[59] +
                   li_Real_s(      -0.085416667304911042535309 ) * liVS[60] +
                   li_Real_s(       0.041875000364308301303851 ) * liVS[61] +
                   li_Real_s(      -0.011793981595470237611067 ) * liVS[62] +
                   li_Real_s(       0.001458333348443252575066 ) * liVS[63];
        liVC[23] = li_Real_s(      -0.000516975309600903321083 ) * liVS[ 0] +
                   li_Real_s(       0.002212301593466886595607 ) * liVS[ 1] +
                   li_Real_s(      -0.004360119065134324123356 ) * liVS[ 2] +
                   li_Real_s(       0.005230379216715458357001 ) * liVS[ 3] +
                   li_Real_s(      -0.004067460344295448226859 ) * liVS[ 4] +
                   li_Real_s(       0.001994047634329310092483 ) * liVS[ 5] +
                   li_Real_s(      -0.000561618170553449674864 ) * liVS[ 6] +
                   li_Real_s(       0.000069444445072439076050 ) * liVS[ 7] +
                   li_Real_s(       0.003618827168504503560831 ) * liVS[ 8] +
                   li_Real_s(      -0.015486111163968413820258 ) * liVS[ 9] +
                   li_Real_s(       0.030520833485494149039940 ) * liVS[10] +
                   li_Real_s(      -0.036612654565599284106625 ) * liVS[11] +
                   li_Real_s(       0.028472222457189780980835 ) * liVS[12] +
                   li_Real_s(      -0.013958333467454462142876 ) * liVS[13] +
                   li_Real_s(       0.003931327202521610678088 ) * liVS[14] +
                   li_Real_s(      -0.000486111116687878985765 ) * liVS[15] +
                   li_Real_s(      -0.010856481507618837012430 ) * liVS[16] +
                   li_Real_s(       0.046458333507627999625633 ) * liVS[17] +
                   li_Real_s(      -0.091562500503982263699854 ) * liVS[18] +
                   li_Real_s(       0.109837963774082197909721 ) * liVS[19] +
                   li_Real_s(      -0.085416667445767480160157 ) * liVS[20] +
                   li_Real_s(       0.041875000444744146788079 ) * liVS[21] +
                   li_Real_s(      -0.011793981620972837642825 ) * liVS[22] +
                   li_Real_s(       0.001458333351887053375151 ) * liVS[23] +
                   li_Real_s(       0.018094135847042280840213 ) * liVS[24] +
                   li_Real_s(      -0.077430555854188884423728 ) * liVS[25] +
                   li_Real_s(       0.152604167531378953626131 ) * liVS[26] +
                   li_Real_s(      -0.183063272996258796787572 ) * liVS[27] +
                   li_Real_s(       0.142361112446440701839734 ) * liVS[28] +
                   li_Real_s(      -0.069791667428376419346492 ) * liVS[29] +
                   li_Real_s(       0.019656636041297366468417 ) * liVS[30] +
                   li_Real_s(      -0.002430555587335240380620 ) * liVS[31] +
                   li_Real_s(      -0.018094135846408371248728 ) * liVS[32] +
                   li_Real_s(       0.077430555850563950737175 ) * liVS[33] +
                   li_Real_s(      -0.152604167520933503565672 ) * liVS[34] +
                   li_Real_s(       0.183063272978245761279936 ) * liVS[35] +
                   li_Real_s(      -0.142361112427622837905972 ) * liVS[36] +
                   li_Real_s(       0.069791667416790589828501 ) * liVS[37] +
                   li_Real_s(      -0.019656636037439563502449 ) * liVS[38] +
                   li_Real_s(       0.002430555586803939682738 ) * liVS[39] +
                   li_Real_s(       0.010856481506859139152255 ) * liVS[40] +
                   li_Real_s(      -0.046458333503896567795444 ) * liVS[41] +
                   li_Real_s(       0.091562500493706455717557 ) * liVS[42] +
                   li_Real_s(      -0.109837963755821305111837 ) * liVS[43] +
                   li_Real_s(       0.085416667425800493762544 ) * liVS[44] +
                   li_Real_s(      -0.041875000431993686378362 ) * liVS[45] +
                   li_Real_s(       0.011793981616633836262409 ) * liVS[46] +
                   li_Real_s(      -0.001458333351288337853546 ) * liVS[47] +
                   li_Real_s(      -0.003618827168500579616328 ) * liVS[48] +
                   li_Real_s(       0.015486111164965093989210 ) * liVS[49] +
                   li_Real_s(      -0.030520833489159578111316 ) * liVS[50] +
                   li_Real_s(       0.036612654571018268812033 ) * liVS[51] +
                   li_Real_s(      -0.028472222461363688728042 ) * liVS[52] +
                   li_Real_s(       0.013958333469259946724161 ) * liVS[53] +
                   li_Real_s(      -0.003931327202966239386939 ) * liVS[54] +
                   li_Real_s(       0.000486111116746769378327 ) * liVS[55] +
                   li_Real_s(       0.000516975309711315000882 ) * liVS[56] +
                   li_Real_s(      -0.002212301594497779848314 ) * liVS[57] +
                   li_Real_s(       0.004360119068436118724974 ) * liVS[58] +
                   li_Real_s(      -0.005230379222093170521468 ) * liVS[59] +
                   li_Real_s(       0.004067460349360028926124 ) * liVS[60] +
                   li_Real_s(      -0.001994047637160601717243 ) * liVS[61] +
                   li_Real_s(       0.000561618171437737109808 ) * liVS[62] +
                   li_Real_s(      -0.000069444445193543591355 ) * liVS[63];
        liVC[24] = li_Real_s(      -1.343055556056689248123348 ) * liVS[ 0] +
                   li_Real_s(       7.088888892184503021098863 ) * liVS[ 1] +
                   li_Real_s(     -16.370833342692854728284146 ) * liVS[ 2] +
                   li_Real_s(      21.611111125902908725038287 ) * liVS[ 3] +
                   li_Real_s(     -17.673611125123532161751427 ) * liVS[ 4] +
                   li_Real_s(       8.933333341287026385657555 ) * liVS[ 5] +
                   li_Real_s(      -2.568055558062077636805043 ) * liVS[ 6] +
                   li_Real_s(       0.322222222560710980232557 ) * liVS[ 7] +
                   li_Real_s(       0.000000003138991012376733 ) * liVS[ 8] +
                   li_Real_s(      -0.000000020712597792665772 ) * liVS[ 9] +
                   li_Real_s(       0.000000058921607346734532 ) * liVS[10] +
                   li_Real_s(      -0.000000093167368260107368 ) * liVS[11] +
                   li_Real_s(       0.000000088240435593722296 ) * liVS[12] +
                   li_Real_s(      -0.000000050057723573998301 ) * liVS[13] +
                   li_Real_s(       0.000000015763949399308340 ) * liVS[14] +
                   li_Real_s(      -0.000000002127293725370436 ) * liVS[15] +
                   li_Real_s(      -0.000000008414898056051889 ) * liVS[16] +
                   li_Real_s(       0.000000055719282368565921 ) * liVS[17] +
                   li_Real_s(      -0.000000158791877983752417 ) * liVS[18] +
                   li_Real_s(       0.000000251235534647522083 ) * liVS[19] +
                   li_Real_s(      -0.000000237911199097038370 ) * liVS[20] +
                   li_Real_s(       0.000000134885629924639097 ) * liVS[21] +
                   li_Real_s(      -0.000000042446097165199777 ) * liVS[22] +
                   li_Real_s(       0.000000005723625361315345 ) * liVS[23] +
                   li_Real_s(       0.000000012473127952253676 ) * liVS[24] +
                   li_Real_s(      -0.000000082904432876541290 ) * liVS[25] +
                   li_Real_s(       0.000000236769583036349537 ) * liVS[26] +
                   li_Real_s(      -0.000000374944163602604393 ) * liVS[27] +
                   li_Real_s(       0.000000355081165846438772 ) * liVS[28] +
                   li_Real_s(      -0.000000201231964874786471 ) * liVS[29] +
                   li_Real_s(       0.000000063284411255730557 ) * liVS[30] +
                   li_Real_s(      -0.000000008527726736839830 ) * liVS[31] +
                   li_Real_s(      -0.000000011039374627883425 ) * liVS[32] +
                   li_Real_s(       0.000000073679499728403632 ) * liVS[33] +
                   li_Real_s(      -0.000000210954192781006126 ) * liVS[34] +
                   li_Real_s(       0.000000334480293102756021 ) * liVS[35] +
                   li_Real_s(      -0.000000316874705800960786 ) * liVS[36] +
                   li_Real_s(       0.000000179547110324104599 ) * liVS[37] +
                   li_Real_s(      -0.000000056439866996393455 ) * liVS[38] +
                   li_Real_s(       0.000000007601237050979527 ) * liVS[39] +
                   li_Real_s(       0.000000005844021334187353 ) * liVS[40] +
                   li_Real_s(      -0.000000039178071843867308 ) * liVS[41] +
                   li_Real_s(       0.000000112493118512006659 ) * liVS[42] +
                   li_Real_s(      -0.000000178645122748696302 ) * liVS[43] +
                   li_Real_s(       0.000000169350775128617370 ) * liVS[44] +
                   li_Real_s(      -0.000000095962538930552661 ) * liVS[45] +
                   li_Real_s(       0.000000030157793485986121 ) * liVS[46] +
                   li_Real_s(      -0.000000004059974937681220 ) * liVS[47] +
                   li_Real_s(      -0.000000001718337347485017 ) * liVS[48] +
                   li_Real_s(       0.000000011572370317591359 ) * liVS[49] +
                   li_Real_s(      -0.000000033329675503841991 ) * liVS[50] +
                   li_Real_s(       0.000000053023814613456994 ) * liVS[51] +
                   li_Real_s(      -0.000000050307509809731676 ) * liVS[52] +
                   li_Real_s(       0.000000028513302194829617 ) * liVS[53] +
                   li_Real_s(      -0.000000008959837800292412 ) * liVS[54] +
                   li_Real_s(       0.000000001205873335473136 ) * liVS[55] +
                   li_Real_s(       0.000000000217301275073673 ) * liVS[56] +
                   li_Real_s(      -0.000000001469942621926555 ) * liVS[57] +
                   li_Real_s(       0.000000004246453307072849 ) * liVS[58] +
                   li_Real_s(      -0.000000006768000696077222 ) * liVS[59] +
                   li_Real_s(       0.000000006427186163418472 ) * liVS[60] +
                   li_Real_s(      -0.000000003643952005221878 ) * liVS[61] +
                   li_Real_s(       0.000000001145024891551486 ) * liVS[62] +
                   li_Real_s(      -0.000000000154070313890881 ) * liVS[63];
        liVC[25] = li_Real_s(       3.482351190452959599497262 ) * liVS[ 0] +
                   li_Real_s(     -18.380476190525556035026966 ) * liVS[ 1] +
                   li_Real_s(      42.447232144231286099511635 ) * liVS[ 2] +
                   li_Real_s(     -56.034523813936630176613107 ) * liVS[ 3] +
                   li_Real_s(      45.825148815656689293973614 ) * liVS[ 4] +
                   li_Real_s(     -23.162857147211980191059411 ) * liVS[ 5] +
                   li_Real_s(       6.658601192030221938011891 ) * liVS[ 6] +
                   li_Real_s(      -0.835476190696610387931287 ) * liVS[ 7] +
                   li_Real_s(      -9.401388891183017904040753 ) * liVS[ 8] +
                   li_Real_s(      49.622222239595686232860317 ) * liVS[ 9] +
                   li_Real_s(    -114.595833391945291168667609 ) * liVS[10] +
                   li_Real_s(     151.277777885595952511721407 ) * liVS[11] +
                   li_Real_s(    -123.715277892482561128417728 ) * liVS[12] +
                   li_Real_s(      62.533333403816250495310669 ) * liVS[13] +
                   li_Real_s(     -17.976388912142695630791422 ) * liVS[14] +
                   li_Real_s(       2.255555558745669486597762 ) * liVS[15] +
                   li_Real_s(      14.102083344949832621750829 ) * liVS[16] +
                   li_Real_s(     -74.433333417940730214468203 ) * liVS[17] +
                   li_Real_s(     171.893750268186835228334530 ) * liVS[18] +
                   li_Real_s(    -226.916667133282828672236064 ) * liVS[19] +
                   li_Real_s(     185.572917142718665672873612 ) * liVS[20] +
                   li_Real_s(     -93.800000284160759633778071 ) * liVS[21] +
                   li_Real_s(      26.964583425376204672829772 ) * liVS[22] +
                   li_Real_s(      -3.383333345847248097015836 ) * liVS[23] +
                   li_Real_s(     -15.668981505026948752856697 ) * liVS[24] +
                   li_Real_s(      82.703703873368169752211543 ) * liVS[25] +
                   li_Real_s(    -190.993056082749063762094011 ) * liVS[26] +
                   li_Real_s(     252.129630529068293753880425 ) * liVS[27] +
                   li_Real_s(    -206.192130532699223977033398 ) * liVS[28] +
                   li_Real_s(     104.222222754918362852549762 ) * liVS[29] +
                   li_Real_s(     -29.960648319305903442000272 ) * liVS[30] +
                   li_Real_s(       3.759259282426073767169328 ) * liVS[31] +
                   li_Real_s(      11.751736136396516485547181 ) * liVS[32] +
                   li_Real_s(     -62.027777959321632295086602 ) * liVS[33] +
                   li_Real_s(     143.244792225743623248490621 ) * liVS[34] +
                   li_Real_s(    -189.097223167015386025013868 ) * liVS[35] +
                   li_Real_s(     154.644098163066445295044105 ) * liVS[36] +
                   li_Real_s(     -78.166667218070656986128597 ) * liVS[37] +
                   li_Real_s(      22.470486287440500916545716 ) * liVS[38] +
                   li_Real_s(      -2.819444468239424850253272 ) * liVS[39] +
                   li_Real_s(      -5.640833348658290447019681 ) * liVS[40] +
                   li_Real_s(      29.773333443271884135583605 ) * liVS[41] +
                   li_Real_s(     -68.757500337011222768524021 ) * liVS[42] +
                   li_Real_s(      90.766667233125005509464245 ) * liVS[43] +
                   li_Real_s(     -74.229167227973903209203854 ) * liVS[44] +
                   li_Real_s(      37.520000327606723544704437 ) * liVS[45] +
                   li_Real_s(     -10.785833437750579122393901 ) * liVS[46] +
                   li_Real_s(       1.353333347390402785492824 ) * liVS[47] +
                   li_Real_s(       1.566898153128924775501218 ) * liVS[48] +
                   li_Real_s(      -8.270370406124094131428137 ) * liVS[49] +
                   li_Real_s(      19.099305664900754209156730 ) * liVS[50] +
                   li_Real_s(     -25.212963146163261285437329 ) * liVS[51] +
                   li_Real_s(      20.619213143918727837444749 ) * liVS[52] +
                   li_Real_s(     -10.422222327533589236736589 ) * liVS[53] +
                   li_Real_s(       2.996064848296867921817466 ) * liVS[54] +
                   li_Real_s(      -0.375925930424315879463393 ) * liVS[55] +
                   li_Real_s(      -0.191865080042276758831576 ) * liVS[56] +
                   li_Real_s(       1.012698417567577280351543 ) * liVS[57] +
                   li_Real_s(      -2.338690491067783483458697 ) * liVS[58] +
                   li_Real_s(       3.087301612181676091495319 ) * liVS[59] +
                   li_Real_s(      -2.524801611826156033657753 ) * liVS[60] +
                   li_Real_s(       1.276190490433900759370545 ) * liVS[61] +
                   li_Real_s(      -0.366865083884707843253636 ) * liVS[62] +
                   li_Real_s(       0.046031746637865467164374 ) * liVS[63];
        liVC[26] = li_Real_s(      -3.499405862170874570438173 ) * liVS[ 0] +
                   li_Real_s(      18.470493814511030450375983 ) * liVS[ 1] +
                   li_Real_s(     -42.655115708212690606160322 ) * liVS[ 2] +
                   li_Real_s(      56.308950571786908767535351 ) * liVS[ 3] +
                   li_Real_s(     -46.049575579241633249694132 ) * liVS[ 4] +
                   li_Real_s(      23.276296276944918872686685 ) * liVS[ 5] +
                   li_Real_s(      -6.691211414120138201155896 ) * liVS[ 6] +
                   li_Real_s(       0.839567900502071751134281 ) * liVS[ 7] +
                   li_Real_s(      14.975069437185766219045036 ) * liVS[ 8] +
                   li_Real_s(     -79.041111069811691436370893 ) * liVS[ 9] +
                   li_Real_s(     182.534791574084493959162501 ) * liVS[10] +
                   li_Real_s(    -240.963888782366836949222488 ) * liVS[11] +
                   li_Real_s(     197.060763820720779904149822 ) * liVS[12] +
                   li_Real_s(     -99.606666641992546828987543 ) * liVS[13] +
                   li_Real_s(      28.633819439424698316543072 ) * liVS[14] +
                   li_Real_s(      -3.592777777244748449447798 ) * liVS[15] +
                   li_Real_s(     -29.513645825789552645801450 ) * liVS[16] +
                   li_Real_s(     155.778333305336389003059594 ) * liVS[17] +
                   li_Real_s(    -359.749062493387896211061161 ) * liVS[18] +
                   li_Real_s(     474.904166767963715756195597 ) * liVS[19] +
                   li_Real_s(    -388.377604351926720482879318 ) * liVS[20] +
                   li_Real_s(     196.310000142207115914061433 ) * liVS[21] +
                   li_Real_s(     -56.433020885006769162828277 ) * liVS[22] +
                   li_Real_s(       7.080833340603817305236589 ) * liVS[23] +
                   li_Real_s(      35.404436730597865334857488 ) * liVS[24] +
                   li_Real_s(    -186.870987706193176336455508 ) * liVS[25] +
                   li_Real_s(     431.553356745206144751136890 ) * liVS[26] +
                   li_Real_s(    -569.692901823633064850582741 ) * liVS[27] +
                   li_Real_s(     465.896026927106277071288787 ) * liVS[28] +
                   li_Real_s(    -235.492593039897258222481469 ) * liVS[29] +
                   li_Real_s(      67.696797989878177759237587 ) * liVS[30] +
                   li_Real_s(      -8.494135823065056101199843 ) * liVS[31] +
                   li_Real_s(     -27.532638899876644700270845 ) * liVS[32] +
                   li_Real_s(     145.322222331746644385930267 ) * liVS[33] +
                   li_Real_s(    -335.602083754128898362978362 ) * liVS[34] +
                   li_Real_s(     443.027778600850922430254286 ) * liVS[35] +
                   li_Real_s(    -362.309028678426841452164808 ) * liVS[36] +
                   li_Real_s(     183.133333891921211034059525 ) * liVS[37] +
                   li_Real_s(     -52.645139072568397864415601 ) * liVS[38] +
                   li_Real_s(       6.605555580482047162149684 ) * liVS[39] +
                   li_Real_s(      13.497708342899898070754716 ) * liVS[40] +
                   li_Real_s(     -71.243333418036527859840135 ) * liVS[41] +
                   li_Real_s(     164.526875302070891393668717 ) * liVS[42] +
                   li_Real_s(    -217.191667230663171039850567 ) * liVS[43] +
                   li_Real_s(     177.619792266003003078367328 ) * liVS[44] +
                   li_Real_s(     -89.780000364939240853345837 ) * liVS[45] +
                   li_Real_s(      25.808958452015104967358639 ) * liVS[46] +
                   li_Real_s(      -3.238333349350009271461204 ) * liVS[47] +
                   li_Real_s(      -3.801593367872357021042262 ) * liVS[48] +
                   li_Real_s(      20.065493858316397535190845 ) * liVS[49] +
                   li_Real_s(     -46.338553348303321399725974 ) * liVS[50] +
                   li_Real_s(      61.171450813623494013882009 ) * liVS[51] +
                   li_Real_s(     -50.026138322734254870738368 ) * liVS[52] +
                   li_Real_s(      25.286296420101095350219111 ) * liVS[53] +
                   li_Real_s(      -7.269023959740025020437315 ) * liVS[54] +
                   li_Real_s(       0.912067906608939438228845 ) * liVS[55] +
                   li_Real_s(       0.470069444988041595934192 ) * liVS[56] +
                   li_Real_s(      -2.481111115637820496715449 ) * liVS[57] +
                   li_Real_s(       5.729791682053779311445396 ) * liVS[58] +
                   li_Real_s(      -7.563888916646163806944969 ) * liVS[59] +
                   li_Real_s(       6.185763917684305113198207 ) * liVS[60] +
                   li_Real_s(      -3.126666683909483879233449 ) * liVS[61] +
                   li_Real_s(       0.898819449987612983932195 ) * liVS[62] +
                   li_Real_s(      -0.112777778520182891952572 ) * liVS[63];
        liVC[27] = li_Real_s(       1.803798223728506400220795 ) * liVS[ 0] +
                   li_Real_s(      -9.520771595168312728674209 ) * liVS[ 1] +
                   li_Real_s(      21.986938632505385271542764 ) * liVS[ 2] +
                   li_Real_s(     -29.024922804997061120957369 ) * liVS[ 3] +
                   li_Real_s(      23.736641560901489356183447 ) * liVS[ 4] +
                   li_Real_s(     -11.997962948497445268003503 ) * liVS[ 5] +
                   li_Real_s(       3.449041276655471222056804 ) * liVS[ 6] +
                   li_Real_s(      -0.432762345128338665745105 ) * liVS[ 7] +
                   li_Real_s(      -9.520771599456068656763819 ) * liVS[ 8] +
                   li_Real_s(      50.252345648714047854355158 ) * liVS[ 9] +
                   li_Real_s(    -116.051018453373018246566062 ) * liVS[10] +
                   li_Real_s(     153.198765362160429504001513 ) * liVS[11] +
                   li_Real_s(    -125.286265392636209980992135 ) * liVS[12] +
                   li_Real_s(      63.327407396234903558251972 ) * liVS[13] +
                   li_Real_s(     -18.204660492405139393667923 ) * liVS[14] +
                   li_Real_s(       2.284197530761019834244507 ) * liVS[15] +
                   li_Real_s(      21.986938652342388422766817 ) * liVS[16] +
                   li_Real_s(    -116.051018504374496842501685 ) * liVS[17] +
                   li_Real_s(     268.004184049581738236156525 ) * liVS[18] +
                   li_Real_s(    -353.791898276010726931417594 ) * liVS[19] +
                   li_Real_s(     289.331742092978174696327187 ) * liVS[20] +
                   li_Real_s(    -146.246111251653019280638546 ) * liVS[21] +
                   li_Real_s(      42.041209540288448920364317 ) * liVS[22] +
                   li_Real_s(      -5.275046303152379323364585 ) * liVS[23] +
                   li_Real_s(     -29.024922842803164257929893 ) * liVS[24] +
                   li_Real_s(     153.198765486668605717568425 ) * liVS[25] +
                   li_Real_s(    -353.791898402167475978785660 ) * liVS[26] +
                   li_Real_s(     467.040124002705056227569003 ) * liVS[27] +
                   li_Real_s(    -381.946374085452248436922673 ) * liVS[28] +
                   li_Real_s(     193.059259660127395363815594 ) * liVS[29] +
                   li_Real_s(     -55.498534084152012724189262 ) * liVS[30] +
                   li_Real_s(       6.963580265073269437436920 ) * liVS[31] +
                   li_Real_s(      23.736641599712584138615057 ) * liVS[32] +
                   li_Real_s(    -125.286265531755120150592120 ) * liVS[33] +
                   li_Real_s(     289.331742274391615410422673 ) * liVS[34] +
                   li_Real_s(    -381.946374183987359174352605 ) * liVS[35] +
                   li_Real_s(     312.356530496059860979585210 ) * liVS[36] +
                   li_Real_s(    -157.884259745758811277482891 ) * liVS[37] +
                   li_Real_s(      45.386815359735976471711183 ) * liVS[38] +
                   li_Real_s(      -5.694830268398646921923500 ) * liVS[39] +
                   li_Real_s(     -11.997962971378996144267148 ) * liVS[40] +
                   li_Real_s(      63.327407481653267495858017 ) * liVS[41] +
                   li_Real_s(    -146.246111374083454848005204 ) * liVS[42] +
                   li_Real_s(     193.059259747432776066489168 ) * liVS[43] +
                   li_Real_s(    -157.884259775682352255898877 ) * liVS[44] +
                   li_Real_s(      79.804444757651339159565396 ) * liVS[45] +
                   li_Real_s(     -22.941296397743258239643183 ) * liVS[46] +
                   li_Real_s(       2.878518532150692976756545 ) * liVS[47] +
                   li_Real_s(       3.449041284023564912786242 ) * liVS[48] +
                   li_Real_s(     -18.204660520635314924220438 ) * liVS[49] +
                   li_Real_s(      42.041209583000579641520744 ) * liVS[50] +
                   li_Real_s(     -55.498534118474054821490427 ) * liVS[51] +
                   li_Real_s(      45.386815375753371881728526 ) * liVS[52] +
                   li_Real_s(     -22.941296401521597658756946 ) * liVS[53] +
                   li_Real_s(       6.594909370288462469034130 ) * liVS[54] +
                   li_Real_s(      -0.827484572435025711456547 ) * liVS[55] +
                   li_Real_s(      -0.432762346139725195826031 ) * liVS[56] +
                   li_Real_s(       2.284197534713124255745242 ) * liVS[57] +
                   li_Real_s(      -5.275046309362103613693762 ) * liVS[58] +
                   li_Real_s(       6.963580270437603303435026 ) * liVS[59] +
                   li_Real_s(      -5.694830271268358501401963 ) * liVS[60] +
                   li_Real_s(       2.878518533067062179497952 ) * liVS[61] +
                   li_Real_s(      -0.827484572563328413252748 ) * liVS[62] +
                   li_Real_s(       0.103827161115830790549808 ) * liVS[63];
        liVC[28] = li_Real_s(      -0.522299382375592813332332 ) * liVS[ 0] +
                   li_Real_s(       2.756790121449199659764417 ) * liVS[ 1] +
                   li_Real_s(      -6.366435180521772707606942 ) * liVS[ 2] +
                   li_Real_s(       8.404320982043202548084082 ) * liVS[ 3] +
                   li_Real_s(      -6.873070983807750167216000 ) * liVS[ 4] +
                   li_Real_s(       3.474074072509296229327447 ) * liVS[ 5] +
                   li_Real_s(      -0.998688271224195567299375 ) * liVS[ 6] +
                   li_Real_s(       0.125308641927450281627898 ) * liVS[ 7] +
                   li_Real_s(       3.105815971831912492007177 ) * liVS[ 8] +
                   li_Real_s(     -16.393055555244767873546152 ) * liVS[ 9] +
                   li_Real_s(      37.857552090313191683890182 ) * liVS[10] +
                   li_Real_s(     -49.975694467728274617002171 ) * liVS[11] +
                   li_Real_s(      40.870225726314075131995196 ) * liVS[12] +
                   li_Real_s(     -20.658333355412608511869621 ) * liVS[13] +
                   li_Real_s(       5.938628479861474573908708 ) * liVS[14] +
                   li_Real_s(      -0.745138889935006432096998 ) * liVS[15] +
                   li_Real_s(      -7.946412039829198192819604 ) * liVS[16] +
                   li_Real_s(      41.942592618581656438436767 ) * liVS[17] +
                   li_Real_s(     -96.860763986355664201255422 ) * liVS[18] +
                   li_Real_s(     127.865740929858276331287925 ) * liVS[19] +
                   li_Real_s(    -104.568865947168802676969790 ) * liVS[20] +
                   li_Real_s(      52.855555683457417615045415 ) * liVS[21] +
                   li_Real_s(     -15.194328745708212125009595 ) * liVS[22] +
                   li_Real_s(       1.906481487164541022139019 ) * liVS[23] +
                   li_Real_s(      11.369338357297209540774929 ) * liVS[24] +
                   li_Real_s(     -60.009413649496615050793480 ) * liVS[25] +
                   li_Real_s(     138.583651855081910753142438 ) * liVS[26] +
                   li_Real_s(    -182.944059067451860300934641 ) * liVS[27] +
                   li_Real_s(     149.612027836459276386449346 ) * liVS[28] +
                   li_Real_s(     -75.623148416225802748158458 ) * liVS[29] +
                   li_Real_s(      21.739303713303513632126851 ) * liVS[30] +
                   li_Real_s(      -2.727700628967451024209367 ) * liVS[31] +
                   li_Real_s(      -9.849074084674271034600679 ) * liVS[32] +
                   li_Real_s(      51.985185267893093907787261 ) * liVS[33] +
                   li_Real_s(    -120.052778048524615428505058 ) * liVS[34] +
                   li_Real_s(     158.481481959461632413876941 ) * liVS[35] +
                   li_Real_s(    -129.606481971546656950522447 ) * liVS[36] +
                   li_Real_s(      65.511111403004150588458288 ) * liVS[37] +
                   li_Real_s(     -18.832407501157639728717186 ) * liVS[38] +
                   li_Real_s(       2.362962975544306232222880 ) * liVS[39] +
                   li_Real_s(       5.176359960616792932341923 ) * liVS[40] +
                   li_Real_s(     -27.321759312257881902041845 ) * liVS[41] +
                   li_Real_s(      63.095920309436081652165740 ) * liVS[42] +
                   li_Real_s(     -83.292824370988554960604233 ) * liVS[43] +
                   li_Real_s(      68.117043125320236640618532 ) * liVS[44] +
                   li_Real_s(     -34.430555733639224058606487 ) * liVS[45] +
                   li_Real_s(       9.897714177273947200319526 ) * liVS[46] +
                   li_Real_s(      -1.241898155761399280549995 ) * liVS[47] +
                   li_Real_s(      -1.529591051720807115543721 ) * liVS[48] +
                   li_Real_s(       8.073456807900230458585611 ) * liVS[49] +
                   li_Real_s(     -18.644560241862365046472405 ) * liVS[50] +
                   li_Real_s(      24.612654418864607919203991 ) * liVS[51] +
                   li_Real_s(     -20.128279419653694048975012 ) * liVS[52] +
                   li_Real_s(      10.174074132110632717740373 ) * liVS[53] +
                   li_Real_s(      -2.924729956745458991917985 ) * liVS[54] +
                   li_Real_s(       0.366975311106855883735989 ) * liVS[55] +
                   li_Real_s(       0.195862268842546427549678 ) * liVS[56] +
                   li_Real_s(      -1.033796298753381748269931 ) * liVS[57] +
                   li_Real_s(       2.387413202241397414127277 ) * liVS[58] +
                   li_Real_s(      -3.151620383773490630119340 ) * liVS[59] +
                   li_Real_s(       2.577401633828038995943643 ) * liVS[60] +
                   li_Real_s(      -1.302777785666856980029138 ) * liVS[61] +
                   li_Real_s(       0.374508104355582460698315 ) * liVS[62] +
                   li_Real_s(      -0.046990741073962283280707 ) * liVS[63];
        liVC[29] = li_Real_s(       0.085806327172591068119800 ) * liVS[ 0] +
                   li_Real_s(      -0.452901234686761355874296 ) * liVS[ 1] +
                   li_Real_s(       1.045914352369064381775843 ) * liVS[ 2] +
                   li_Real_s(      -1.380709877668742535661295 ) * liVS[ 3] +
                   li_Real_s(       1.129147377865502122062935 ) * liVS[ 4] +
                   li_Real_s(      -0.570740741595010092623852 ) * liVS[ 5] +
                   li_Real_s(       0.164070216334815555114801 ) * liVS[ 6] +
                   li_Real_s(      -0.020586419791481347374429 ) * liVS[ 7] +
                   li_Real_s(      -0.550279707230885151147959 ) * liVS[ 8] +
                   li_Real_s(       2.904475311962873007587405 ) * liVS[ 9] +
                   li_Real_s(      -6.707494223828326163072688 ) * liVS[10] +
                   li_Real_s(       8.854552488528405262968590 ) * liVS[11] +
                   li_Real_s(      -7.241271239263936898566953 ) * liVS[12] +
                   li_Real_s(       3.660185197295365622238705 ) * liVS[13] +
                   li_Real_s(      -1.052189432928335666161956 ) * liVS[14] +
                   li_Real_s(       0.132021605464842650690116 ) * liVS[15] +
                   li_Real_s(       1.510937501997897669525628 ) * liVS[16] +
                   li_Real_s(      -7.975000014615342358581529 ) * liVS[17] +
                   li_Real_s(      18.417187545925614244879398 ) * liVS[18] +
                   li_Real_s(     -24.312500079047353551686683 ) * liVS[19] +
                   li_Real_s(      19.882812579830773103140018 ) * liVS[20] +
                   li_Real_s(     -10.050000047171595340955719 ) * liVS[21] +
                   li_Real_s(       2.889062515108548723219428 ) * liVS[22] +
                   li_Real_s(      -0.362500002028545154075800 ) * liVS[23] +
                   li_Real_s(      -2.303713352690753168872106 ) * liVS[24] +
                   li_Real_s(      12.159413608637965609204912 ) * liVS[25] +
                   li_Real_s(     -28.080526708175277406098758 ) * liVS[26] +
                   li_Real_s(      37.069058790859997998268227 ) * liVS[27] +
                   li_Real_s(     -30.315152540556919547043435 ) * liVS[28] +
                   li_Real_s(      15.323148235209792744626611 ) * liVS[29] +
                   li_Real_s(      -4.404928654287225242569548 ) * liVS[30] +
                   li_Real_s(       0.552700621002449210550367 ) * liVS[31] +
                   li_Real_s(       2.107851084418637555017995 ) * liVS[32] +
                   li_Real_s(     -11.125617313942832709017239 ) * liVS[33] +
                   li_Real_s(      25.693113517882224527966173 ) * liVS[34] +
                   li_Real_s(     -33.917438426225686498582945 ) * liVS[35] +
                   li_Real_s(      27.737750924842440980455649 ) * liVS[36] +
                   li_Real_s(     -14.020370459689248221479829 ) * liVS[37] +
                   li_Real_s(       4.030420553051058618621028 ) * liVS[38] +
                   li_Real_s(      -0.505709880336597805694510 ) * liVS[39] +
                   li_Real_s(      -1.158385419186370945965336 ) * liVS[40] +
                   li_Real_s(       6.114166684720291122800973 ) * liVS[41] +
                   li_Real_s(     -14.119843805055488417110610 ) * liVS[42] +
                   li_Real_s(      18.639583425406936356694132 ) * liVS[43] +
                   li_Real_s(     -15.243489674156489144252191 ) * liVS[44] +
                   li_Real_s(       7.705000052737581484052498 ) * liVS[45] +
                   li_Real_s(      -2.214947933361505860716534 ) * liVS[46] +
                   li_Real_s(       0.277916668895055229970836 ) * liVS[47] +
                   li_Real_s(       0.354417439088898333920952 ) * liVS[48] +
                   li_Real_s(      -1.870679018192556553046302 ) * liVS[49] +
                   li_Real_s(       4.320081036290927656295935 ) * liVS[50] +
                   li_Real_s(      -5.702932128382876442174165 ) * liVS[51] +
                   li_Real_s(       4.663869627885993907057127 ) * liVS[52] +
                   li_Real_s(      -2.357407424268489926078018 ) * liVS[53] +
                   li_Real_s(       0.677681332484855403208712 ) * liVS[54] +
                   li_Real_s(      -0.085030864906753933496475 ) * liVS[55] +
                   li_Real_s(      -0.046633873567724748454566 ) * liVS[56] +
                   li_Real_s(       0.246141976101841936097614 ) * liVS[57] +
                   li_Real_s(      -0.568431715369748014055062 ) * liVS[58] +
                   li_Real_s(       0.750385806471225436098393 ) * liVS[59] +
                   li_Real_s(      -0.613667056395392318535187 ) * liVS[60] +
                   li_Real_s(       0.310185187453656252554168 ) * liVS[61] +
                   li_Real_s(      -0.089168596393842891600912 ) * liVS[62] +
                   li_Real_s(       0.011188271699959173588468 ) * liVS[63];
        liVC[30] = li_Real_s(      -0.007461419763727938914144 ) * liVS[ 0] +
                   li_Real_s(       0.039382716122031258265679 ) * liVS[ 1] +
                   li_Real_s(      -0.090949074294283782649018 ) * liVS[ 2] +
                   li_Real_s(       0.120061728767888453717205 ) * liVS[ 3] +
                   li_Real_s(      -0.098186728768878683837329 ) * liVS[ 4] +
                   li_Real_s(       0.049629629849676115682655 ) * liVS[ 5] +
                   li_Real_s(      -0.014266975378864454881978 ) * liVS[ 6] +
                   li_Real_s(       0.001790123466158866083475 ) * liVS[ 7] +
                   li_Real_s(       0.050364583438886501021159 ) * liVS[ 8] +
                   li_Real_s(      -0.265833334063509507672052 ) * liVS[ 9] +
                   li_Real_s(       0.613906252194892698703654 ) * liVS[10] +
                   li_Real_s(      -0.810416670323796495267743 ) * liVS[11] +
                   li_Real_s(       0.662760420278308171049275 ) * liVS[12] +
                   li_Real_s(      -0.335000002103797267505314 ) * liVS[13] +
                   li_Real_s(       0.096302084001789634015722 ) * liVS[14] +
                   li_Real_s(      -0.012083333422773900878155 ) * liVS[15] +
                   li_Real_s(      -0.145497685556816636065491 ) * liVS[16] +
                   li_Real_s(       0.767962965540332387348599 ) * liVS[17] +
                   li_Real_s(      -1.773506952150949445012884 ) * liVS[18] +
                   li_Real_s(       2.341203716443913229738882 ) * liVS[19] +
                   li_Real_s(      -1.914641216193787132127113 ) * liVS[20] +
                   li_Real_s(       0.967777785014333646884666 ) * liVS[21] +
                   li_Real_s(      -0.278206020811220700039712 ) * liVS[22] +
                   li_Real_s(       0.034907407714194427228449 ) * liVS[23] +
                   li_Real_s(       0.233169367943109495655563 ) * liVS[24] +
                   li_Real_s(      -1.230709881115719772637362 ) * liVS[25] +
                   li_Real_s(       2.842158578437259386362257 ) * liVS[26] +
                   li_Real_s(      -3.751929034758325620657615 ) * liVS[27] +
                   li_Real_s(       3.068335284220343517347374 ) * liVS[28] +
                   li_Real_s(      -1.550925938555261396345486 ) * liVS[29] +
                   li_Real_s(       0.445842982386966735219858 ) * liVS[30] +
                   li_Real_s(      -0.055941358558369562448132 ) * liVS[31] +
                   li_Real_s(      -0.223842593257656652383503 ) * liVS[32] +
                   li_Real_s(       1.181481486092668786724857 ) * liVS[33] +
                   li_Real_s(      -2.728472235920799793262859 ) * liVS[34] +
                   li_Real_s(       3.601851874310947554391760 ) * liVS[35] +
                   li_Real_s(      -2.945601873699510431237059 ) * liVS[36] +
                   li_Real_s(       1.488888901467460001981635 ) * liVS[37] +
                   li_Real_s(      -0.428009263226433289695194 ) * liVS[38] +
                   li_Real_s(       0.053703704233323268368849 ) * liVS[39] +
                   li_Real_s(       0.128709491130769793088007 ) * liVS[40] +
                   li_Real_s(      -0.679351854554706369171413 ) * liVS[41] +
                   li_Real_s(       1.568871535789457150045223 ) * liVS[42] +
                   li_Real_s(      -2.071064827914754324922342 ) * liVS[43] +
                   li_Real_s(       1.693721077524833473759713 ) * liVS[44] +
                   li_Real_s(      -0.856111118411553295715066 ) * liVS[45] +
                   li_Real_s(       0.246105326371831023557490 ) * liVS[46] +
                   li_Real_s(      -0.030879629935877596358385 ) * liVS[47] +
                   li_Real_s(      -0.041037808766341665744903 ) * liVS[48] +
                   li_Real_s(       0.216604939133362955905326 ) * liVS[49] +
                   li_Real_s(      -0.500219909957950292778150 ) * liVS[50] +
                   li_Real_s(       0.660339510335091306458821 ) * liVS[51] +
                   li_Real_s(      -0.540027010203131485610584 ) * liVS[52] +
                   li_Real_s(       0.272962965273430113732900 ) * liVS[53] +
                   li_Real_s(      -0.078468364923388822518291 ) * liVS[54] +
                   li_Real_s(       0.009845679108928084843910 ) * liVS[55] +
                   li_Real_s(       0.005596064831542513218210 ) * liVS[56] +
                   li_Real_s(      -0.029537037152982850707339 ) * liVS[57] +
                   li_Real_s(       0.068211805898404531678381 ) * liVS[58] +
                   li_Real_s(      -0.090046296855038843176544 ) * liVS[59] +
                   li_Real_s(       0.073640046836519812423205 ) * liVS[60] +
                   li_Real_s(      -0.037222222531434083492297 ) * liVS[61] +
                   li_Real_s(       0.010700231578464308723753 ) * liVS[62] +
                   li_Real_s(      -0.001342592605476422562560 ) * liVS[63];
        liVC[31] = li_Real_s(       0.000266479277714073981009 ) * liVS[ 0] +
                   li_Real_s(      -0.001406525578613433691277 ) * liVS[ 1] +
                   li_Real_s(       0.003248181232681234681392 ) * liVS[ 2] +
                   li_Real_s(      -0.004287918896845113736482 ) * liVS[ 3] +
                   li_Real_s(       0.003506668896086813369539 ) * liVS[ 4] +
                   li_Real_s(      -0.001772486786783405987578 ) * liVS[ 5] +
                   li_Real_s(       0.000509534836959788656197 ) * liVS[ 6] +
                   li_Real_s(      -0.000063932981199884414414 ) * liVS[ 7] +
                   li_Real_s(      -0.001865354945057830227118 ) * liVS[ 8] +
                   li_Real_s(       0.009845679058053537213713 ) * liVS[ 9] +
                   li_Real_s(      -0.022737268652008667257913 ) * liVS[10] +
                   li_Real_s(       0.030015432315555659603667 ) * liVS[11] +
                   li_Real_s(      -0.024546682308623040541518 ) * liVS[12] +
                   li_Real_s(       0.012407407527987922984947 ) * liVS[13] +
                   li_Real_s(      -0.003566743865181468825520 ) * liVS[14] +
                   li_Real_s(       0.000447530869273890519189 ) * liVS[15] +
                   li_Real_s(       0.005596064836860162317045 ) * liVS[16] +
                   li_Real_s(      -0.029537037186436188029859 ) * liVS[17] +
                   li_Real_s(       0.068211805992281798749488 ) * liVS[18] +
                   li_Real_s(      -0.090046297004434991406718 ) * liVS[19] +
                   li_Real_s(       0.073640046980263565434832 ) * liVS[20] +
                   li_Real_s(      -0.037222222614485128699968 ) * liVS[21] +
                   li_Real_s(       0.010700231605048862726370 ) * liVS[22] +
                   li_Real_s(      -0.001342592609098081091190 ) * liVS[23] +
                   li_Real_s(      -0.009326774728885350862484 ) * liVS[24] +
                   li_Real_s(       0.049228395316659449842955 ) * liVS[25] +
                   li_Real_s(      -0.113686343337457751445640 ) * liVS[26] +
                   li_Real_s(       0.150077161699519290305460 ) * liVS[27] +
                   li_Real_s(      -0.122733411656039625992065 ) * liVS[28] +
                   li_Real_s(       0.062037037702392952565056 ) * liVS[29] +
                   li_Real_s(      -0.017833719345119833366020 ) * liVS[30] +
                   li_Real_s(       0.002237654348930764435649 ) * liVS[31] +
                   li_Real_s(       0.009326774728356954091701 ) * liVS[32] +
                   li_Real_s(      -0.049228395313263062504916 ) * liVS[33] +
                   li_Real_s(       0.113686343326748914850377 ) * liVS[34] +
                   li_Real_s(      -0.150077161680059134596377 ) * liVS[35] +
                   li_Real_s(       0.122733411635093839220367 ) * liVS[36] +
                   li_Real_s(      -0.062037037689234700299501 ) * liVS[37] +
                   li_Real_s(       0.017833719340661628727229 ) * liVS[38] +
                   li_Real_s(      -0.002237654348304460305563 ) * liVS[39] +
                   li_Real_s(      -0.005596064836216815829850 ) * liVS[40] +
                   li_Real_s(       0.029537037182563417869741 ) * liVS[41] +
                   li_Real_s(      -0.068211805979811329647688 ) * liVS[42] +
                   li_Real_s(       0.090046296980659981668005 ) * liVS[43] +
                   li_Real_s(      -0.073640046953563700893319 ) * liVS[44] +
                   li_Real_s(       0.037222222597188110715383 ) * liVS[45] +
                   li_Real_s(      -0.010700231599070277044294 ) * liVS[46] +
                   li_Real_s(       0.001342592608250633978706 ) * liVS[47] +
                   li_Real_s(       0.001865354945041294842945 ) * liVS[48] +
                   li_Real_s(      -0.009845679058385268384024 ) * liVS[49] +
                   li_Real_s(       0.022737268652615660879945 ) * liVS[50] +
                   li_Real_s(      -0.030015432314787773848686 ) * liVS[51] +
                   li_Real_s(       0.024546682305932335244059 ) * liVS[52] +
                   li_Real_s(      -0.012407407525417434024373 ) * liVS[53] +
                   li_Real_s(       0.003566743864112133233224 ) * liVS[54] +
                   li_Real_s(      -0.000447530869110958351431 ) * liVS[55] +
                   li_Real_s(      -0.000266479277802947334131 ) * liVS[56] +
                   li_Real_s(       0.001406525579362133404615 ) * liVS[57] +
                   li_Real_s(      -0.003248181234890155227868 ) * liVS[58] +
                   li_Real_s(       0.004287918900153675494380 ) * liVS[59] +
                   li_Real_s(      -0.003506668898936606687533 ) * liVS[60] +
                   li_Real_s(       0.001772486788236725223367 ) * liVS[61] +
                   li_Real_s(      -0.000509534837376340865589 ) * liVS[62] +
                   li_Real_s(       0.000063932981253574105995 ) * liVS[63];
        liVC[32] = li_Real_s(       0.388888889100490908745655 ) * liVS[ 0] +
                   li_Real_s(      -2.312500001395954463134785 ) * liVS[ 1] +
                   li_Real_s(       5.916666670640140957004860 ) * liVS[ 2] +
                   li_Real_s(      -8.465277784067298583181582 ) * liVS[ 3] +
                   li_Real_s(       7.333333339298000908001995 ) * liVS[ 4] +
                   li_Real_s(      -3.854166670055136023620435 ) * liVS[ 5] +
                   li_Real_s(       1.138888889957485384130109 ) * liVS[ 6] +
                   li_Real_s(      -0.145833333477730753280355 ) * liVS[ 7] +
                   li_Real_s(      -0.000000001318526462882907 ) * liVS[ 8] +
                   li_Real_s(       0.000000008729077623809223 ) * liVS[ 9] +
                   li_Real_s(      -0.000000024892495330046478 ) * liVS[10] +
                   li_Real_s(       0.000000039428927736295250 ) * liVS[11] +
                   li_Real_s(      -0.000000037390948024019403 ) * liVS[12] +
                   li_Real_s(       0.000000021232072965259153 ) * liVS[13] +
                   li_Real_s(      -0.000000006691810045693249 ) * liVS[14] +
                   li_Real_s(       0.000000000903701537278418 ) * liVS[15] +
                   li_Real_s(       0.000000003521225428377005 ) * liVS[16] +
                   li_Real_s(      -0.000000023395563132729675 ) * liVS[17] +
                   li_Real_s(       0.000000066844978801824721 ) * liVS[18] +
                   li_Real_s(      -0.000000105957727649971371 ) * liVS[19] +
                   li_Real_s(       0.000000100476701048489441 ) * liVS[20] +
                   li_Real_s(      -0.000000057027426387825101 ) * liVS[21] +
                   li_Real_s(       0.000000017961821590804010 ) * liVS[22] +
                   li_Real_s(      -0.000000002424009698969016 ) * liVS[23] +
                   li_Real_s(      -0.000000005206372832019447 ) * liVS[24] +
                   li_Real_s(       0.000000034725791759511268 ) * liVS[25] +
                   li_Real_s(      -0.000000099436467669487739 ) * liVS[26] +
                   li_Real_s(       0.000000157772641704918514 ) * liVS[27] +
                   li_Real_s(      -0.000000149632281419757196 ) * liVS[28] +
                   li_Real_s(       0.000000084897099509945291 ) * liVS[29] +
                   li_Real_s(      -0.000000026724711402014430 ) * liVS[30] +
                   li_Real_s(       0.000000003604300348903736 ) * liVS[31] +
                   li_Real_s(       0.000000004601004629945382 ) * liVS[32] +
                   li_Real_s(      -0.000000030817062157866954 ) * liVS[33] +
                   li_Real_s(       0.000000088470217033108708 ) * liVS[34] +
                   li_Real_s(      -0.000000140554354298416817 ) * liVS[35] +
                   li_Real_s(       0.000000133356465422834783 ) * liVS[36] +
                   li_Real_s(      -0.000000075652366057775891 ) * liVS[37] +
                   li_Real_s(       0.000000023804971749735640 ) * liVS[38] +
                   li_Real_s(      -0.000000003208876321564835 ) * liVS[39] +
                   li_Real_s(      -0.000000002433591058652604 ) * liVS[40] +
                   li_Real_s(       0.000000016373213205211032 ) * liVS[41] +
                   li_Real_s(      -0.000000047140365811608889 ) * liVS[42] +
                   li_Real_s(       0.000000075012582676302806 ) * liVS[43] +
                   li_Real_s(      -0.000000071218893809620157 ) * liVS[44] +
                   li_Real_s(       0.000000040405495573218938 ) * liVS[45] +
                   li_Real_s(      -0.000000012711260026691888 ) * liVS[46] +
                   li_Real_s(       0.000000001712819251840733 ) * liVS[47] +
                   li_Real_s(       0.000000000715197191643736 ) * liVS[48] +
                   li_Real_s(      -0.000000004834097478892070 ) * liVS[49] +
                   li_Real_s(       0.000000013960831282321743 ) * liVS[50] +
                   li_Real_s(      -0.000000022255434326055297 ) * liVS[51] +
                   li_Real_s(       0.000000021148183180982275 ) * liVS[52] +
                   li_Real_s(      -0.000000012001344896499440 ) * liVS[53] +
                   li_Real_s(       0.000000003775244211541149 ) * liVS[54] +
                   li_Real_s(      -0.000000000508579165042096 ) * liVS[55] +
                   li_Real_s(      -0.000000000090411310540812 ) * liVS[56] +
                   li_Real_s(       0.000000000613854956961933 ) * liVS[57] +
                   li_Real_s(      -0.000000001778254430525651 ) * liVS[58] +
                   li_Real_s(       0.000000002840033272980638 ) * liVS[59] +
                   li_Real_s(      -0.000000002701288483223506 ) * liVS[60] +
                   li_Real_s(       0.000000001533478020020317 ) * liVS[61] +
                   li_Real_s(      -0.000000000482387492667074 ) * liVS[62] +
                   li_Real_s(       0.000000000064975466994241 ) * liVS[63];
        liVC[33] = li_Real_s(      -1.008333333217137806059327 ) * liVS[ 0] +
                   li_Real_s(       5.995982142261993885767879 ) * liVS[ 1] +
                   li_Real_s(     -15.341071427597933407582786 ) * liVS[ 2] +
                   li_Real_s(      21.949255952040942929670564 ) * liVS[ 3] +
                   li_Real_s(     -19.014285714985390995934722 ) * liVS[ 4] +
                   li_Real_s(       9.993303572290407998934825 ) * liVS[ 5] +
                   li_Real_s(      -2.952976190851174909113297 ) * liVS[ 6] +
                   li_Real_s(       0.378125000058307847439210 ) * liVS[ 7] +
                   li_Real_s(       2.722222222421663673230796 ) * liVS[ 8] +
                   li_Real_s(     -16.187500002757062134151056 ) * liVS[ 9] +
                   li_Real_s(      41.416666679543283180464641 ) * liVS[10] +
                   li_Real_s(     -59.256944472655995070908830 ) * liVS[11] +
                   li_Real_s(      51.333333366493434368749149 ) * liVS[12] +
                   li_Real_s(     -26.979166688307355315146197 ) * liVS[13] +
                   li_Real_s(       7.972222229636193091550922 ) * liVS[14] +
                   li_Real_s(      -1.020833334374174228287302 ) * liVS[15] +
                   li_Real_s(      -4.083333335952438858384994 ) * liVS[16] +
                   li_Real_s(      24.281250021927618831796281 ) * liVS[17] +
                   li_Real_s(     -62.125000076946960803070397 ) * liVS[18] +
                   li_Real_s(      88.885416810001061094226316 ) * liVS[19] +
                   li_Real_s(     -77.000000152987468027276918 ) * liVS[20] +
                   li_Real_s(      40.468750094110390591595205 ) * liVS[21] +
                   li_Real_s(     -11.958333364435171475292918 ) * liVS[22] +
                   li_Real_s(       1.531250004282981080905301 ) * liVS[23] +
                   li_Real_s(       4.537037043333825181434804 ) * liVS[24] +
                   li_Real_s(     -26.979166716138795578672216 ) * liVS[25] +
                   li_Real_s(      69.027777941601030420315510 ) * liVS[26] +
                   li_Real_s(     -98.761574366200591157394229 ) * liVS[27] +
                   li_Real_s(      85.555555857869677538474207 ) * liVS[28] +
                   li_Real_s(     -44.965277959832178567012306 ) * liVS[29] +
                   li_Real_s(      13.287037096357433796356418 ) * liVS[30] +
                   li_Real_s(      -1.701388896990318810864551 ) * liVS[31] +
                   li_Real_s(      -3.402777785048982650550897 ) * liVS[32] +
                   li_Real_s(      20.234375055860844838662160 ) * liVS[33] +
                   li_Real_s(     -51.770833514000891284467798 ) * liVS[34] +
                   li_Real_s(      74.071180871482127372473769 ) * liVS[35] +
                   li_Real_s(     -64.166666988782083080877783 ) * liVS[36] +
                   li_Real_s(      33.723958525206789715866762 ) * liVS[37] +
                   li_Real_s(      -9.965277839814675076013373 ) * liVS[38] +
                   li_Real_s(       1.276041675096863059479801 ) * liVS[39] +
                   li_Real_s(       1.633333337924263162221905 ) * liVS[40] +
                   li_Real_s(      -9.712500034915954927328130 ) * liVS[41] +
                   li_Real_s(      24.850000111565112348444018 ) * liVS[42] +
                   li_Real_s(     -35.554166859660476518456562 ) * liVS[43] +
                   li_Real_s(      30.800000195071675079816487 ) * liVS[44] +
                   li_Real_s(     -16.187500115418103519004944 ) * liVS[45] +
                   li_Real_s(       4.783333370460738365181896 ) * liVS[46] +
                   li_Real_s(      -0.612500005027247329536522 ) * liVS[47] +
                   li_Real_s(      -0.453703705235593446332132 ) * liVS[48] +
                   li_Real_s(       2.697916678264211665094763 ) * liVS[49] +
                   li_Real_s(      -6.902777814584112547890982 ) * liVS[50] +
                   li_Real_s(       9.876157470663287085699267 ) * liVS[51] +
                   li_Real_s(      -8.555555619138246470356535 ) * liVS[52] +
                   li_Real_s(       4.496527815227041813272990 ) * liVS[53] +
                   li_Real_s(      -1.328703715706119226069859 ) * liVS[54] +
                   li_Real_s(       0.170138890509523577065920 ) * liVS[55] +
                   li_Real_s(       0.055555555767533348898723 ) * liVS[56] +
                   li_Real_s(      -0.330357144459898277588650 ) * liVS[57] +
                   li_Real_s(       0.845238100306019646268396 ) * liVS[58] +
                   li_Real_s(      -1.209325405500845107553687 ) * liVS[59] +
                   li_Real_s(       1.047619056307830476271192 ) * liVS[60] +
                   li_Real_s(      -0.550595243196618122638597 ) * liVS[61] +
                   li_Real_s(       0.162698414328882989821068 ) * liVS[62] +
                   li_Real_s(      -0.020833333552925381582099 ) * liVS[63];
        liVC[34] = li_Real_s(       1.013271603932793141211732 ) * liVS[ 0] +
                   li_Real_s(      -6.025347216059948607380647 ) * liVS[ 1] +
                   li_Real_s(      15.416203688007840355567168 ) * liVS[ 2] +
                   li_Real_s(     -22.056751521370983937231358 ) * liVS[ 3] +
                   li_Real_s(      19.107407389226004568172357 ) * liVS[ 4] +
                   li_Real_s(     -10.042245361181187490728917 ) * liVS[ 5] +
                   li_Real_s(       2.967438268958279223852514 ) * liVS[ 6] +
                   li_Real_s(      -0.379976851512809687960726 ) * liVS[ 7] +
                   li_Real_s(      -4.336111106919069868581573 ) * liVS[ 8] +
                   li_Real_s(      25.784374976074840191131443 ) * liVS[ 9] +
                   li_Real_s(     -65.970833277907829028663400 ) * liVS[10] +
                   li_Real_s(      94.387847154048813536064699 ) * liVS[11] +
                   li_Real_s(     -81.766666617959074869759206 ) * liVS[12] +
                   li_Real_s(      42.973958312627971167785290 ) * liVS[13] +
                   li_Real_s(     -12.698611106046239882516602 ) * liVS[14] +
                   li_Real_s(       1.626041666080595859966706 ) * liVS[15] +
                   li_Real_s(       8.545833326711999688996002 ) * liVS[16] +
                   li_Real_s(     -50.817187467945380774381192 ) * liVS[17] +
                   li_Real_s(     130.018749944980953614503960 ) * liVS[18] +
                   li_Real_s(    -186.024479133349274206921109 ) * liVS[19] +
                   li_Real_s(     161.150000011147739087391528 ) * liVS[20] +
                   li_Real_s(     -84.695312524477643023601559 ) * liVS[21] +
                   li_Real_s(      25.027083344679077470118500 ) * liVS[22] +
                   li_Real_s(      -3.204687501747471856106131 ) * liVS[23] +
                   li_Real_s(     -10.251543205253284440914285 ) * liVS[24] +
                   li_Real_s(      60.960069433262759730496327 ) * liVS[25] +
                   li_Real_s(    -155.969907432265387114966870 ) * liVS[26] +
                   li_Real_s(     223.154128207475793033154332 ) * liVS[27] +
                   li_Real_s(    -193.314814993225070338667138 ) * liVS[28] +
                   li_Real_s(     101.600115868384321515804913 ) * liVS[29] +
                   li_Real_s(     -30.022376588444963374513463 ) * liVS[30] +
                   li_Real_s(       3.844328710065707532805845 ) * liVS[31] +
                   li_Real_s(       7.972222221647854212278617 ) * liVS[32] +
                   li_Real_s(     -47.406250015004275155661162 ) * liVS[33] +
                   li_Real_s(     121.291666762101201015866536 ) * liVS[34] +
                   li_Real_s(    -173.538194669708957462717080 ) * liVS[35] +
                   li_Real_s(     150.333333603661998267853050 ) * liVS[36] +
                   li_Real_s(     -79.010416843069663173082517 ) * liVS[37] +
                   li_Real_s(      23.347222281934453036456034 ) * liVS[38] +
                   li_Real_s(      -2.989583341562635609989229 ) * liVS[39] +
                   li_Real_s(      -3.908333334498493627506832 ) * liVS[40] +
                   li_Real_s(      23.240625018498576537240297 ) * liVS[41] +
                   li_Real_s(     -59.462500081924545725087228 ) * liVS[42] +
                   li_Real_s(      85.076041836697754661145154 ) * liVS[43] +
                   li_Real_s(     -73.700000191554295270179864 ) * liVS[44] +
                   li_Real_s(      38.734375120706800998959807 ) * liVS[45] +
                   li_Real_s(     -11.445833373391181453371246 ) * liVS[46] +
                   li_Real_s(       1.465625005465396535342393 ) * liVS[47] +
                   li_Real_s(       1.100771605622242077515693 ) * liVS[48] +
                   li_Real_s(      -6.545659730162673994868783 ) * liVS[49] +
                   li_Real_s(      16.747453735317151313211070 ) * liVS[50] +
                   li_Real_s(     -23.961439105453123943334504 ) * liVS[51] +
                   li_Real_s(      20.757407475452822609440773 ) * liVS[52] +
                   li_Real_s(     -10.909432912471523025033093 ) * liVS[53] +
                   li_Real_s(       3.223688285419598287262488 ) * liVS[54] +
                   li_Real_s(      -0.412789353724484442409448 ) * liVS[55] +
                   li_Real_s(      -0.136111111229467951488914 ) * liVS[56] +
                   li_Real_s(       0.809375001244054148585860 ) * liVS[57] +
                   li_Real_s(      -2.070833338063037487586371 ) * liVS[58] +
                   li_Real_s(       2.962847231293800120965898 ) * liVS[59] +
                   li_Real_s(      -2.566666676423764670289529 ) * liVS[60] +
                   li_Real_s(       1.348958339306130627122116 ) * liVS[61] +
                   li_Real_s(      -0.398611113056849930558201 ) * liVS[62] +
                   li_Real_s(       0.051041666929123152840475 ) * liVS[63];
        liVC[35] = li_Real_s(      -0.522299381972260334805469 ) * liVS[ 0] +
                   li_Real_s(       3.105815967720090498005447 ) * liVS[ 1] +
                   li_Real_s(      -7.946412025730332118200749 ) * liVS[ 2] +
                   li_Real_s(      11.369338333288169451407157 ) * liVS[ 3] +
                   li_Real_s(      -9.849074061408721547650202 ) * liVS[ 4] +
                   li_Real_s(       5.176359947397898508825165 ) * liVS[ 5] +
                   li_Real_s(      -1.529591047582557372663814 ) * liVS[ 6] +
                   li_Real_s(       0.195862268287829266455446 ) * liVS[ 7] +
                   li_Real_s(       2.756790120520861364639131 ) * liVS[ 8] +
                   li_Real_s(     -16.393055539395913200451105 ) * liVS[ 9] +
                   li_Real_s(      41.942592557065239589064731 ) * liVS[10] +
                   li_Real_s(     -60.009413539913907698064577 ) * liVS[11] +
                   li_Real_s(      51.985185159744048633001512 ) * liVS[12] +
                   li_Real_s(     -27.321759250292242171553880 ) * liVS[13] +
                   li_Real_s(       8.073456788400278583139880 ) * liVS[14] +
                   li_Real_s(      -1.033796296128357994348335 ) * liVS[15] +
                   li_Real_s(      -6.366435181071324223012198 ) * liVS[16] +
                   li_Real_s(      37.857552066191800577144022 ) * liVS[17] +
                   li_Real_s(     -96.860763870186517010552052 ) * liVS[18] +
                   li_Real_s(     138.583651634773673322342802 ) * liVS[19] +
                   li_Real_s(    -120.052777826008338024621480 ) * liVS[20] +
                   li_Real_s(      63.095920180669665455752693 ) * liVS[21] +
                   li_Real_s(     -18.644560201099682927861068 ) * liVS[22] +
                   li_Real_s(       2.387413196730719278093602 ) * liVS[23] +
                   li_Real_s(       8.404320985780628916472779 ) * liVS[24] +
                   li_Real_s(     -49.975694448653811718941142 ) * liVS[25] +
                   li_Real_s(     127.865740799904784807949909 ) * liVS[26] +
                   li_Real_s(    -182.944058802705200150739984 ) * liVS[27] +
                   li_Real_s(     158.481481685532486380907358 ) * liVS[28] +
                   li_Real_s(     -83.292824210885925140246400 ) * liVS[29] +
                   li_Real_s(      24.612654367878398886659852 ) * liVS[30] +
                   li_Real_s(      -3.151620376851469007561946 ) * liVS[31] +
                   li_Real_s(      -6.873070988867851838222123 ) * liVS[32] +
                   li_Real_s(      40.870225718508621071123343 ) * liVS[33] +
                   li_Real_s(    -104.568865853105734231576207 ) * liVS[34] +
                   li_Real_s(     149.612027630789384602394421 ) * liVS[35] +
                   li_Real_s(    -129.606481754012321516711381 ) * liVS[36] +
                   li_Real_s(      68.117042997028164563744213 ) * liVS[37] +
                   li_Real_s(     -20.128279378567452084780598 ) * liVS[38] +
                   li_Real_s(       2.577401628227178775887296 ) * liVS[39] +
                   li_Real_s(       3.474074075919745041574060 ) * liVS[40] +
                   li_Real_s(     -20.658333354712461016333691 ) * liVS[41] +
                   li_Real_s(      52.855555640726102240023465 ) * liVS[42] +
                   li_Real_s(     -75.623148316086471254493517 ) * liVS[43] +
                   li_Real_s(      65.511111294969126106479962 ) * liVS[44] +
                   li_Real_s(     -34.430555669416222031031793 ) * liVS[45] +
                   li_Real_s(      10.174074111438546452745868 ) * liVS[46] +
                   li_Real_s(      -1.302777782838340669968602 ) * liVS[47] +
                   li_Real_s(      -0.998688272426747403187619 ) * liVS[48] +
                   li_Real_s(       5.938628480552589294916288 ) * liVS[49] +
                   li_Real_s(     -15.194328734890717669259175 ) * liVS[50] +
                   li_Real_s(      21.739303685914293851055845 ) * liVS[51] +
                   li_Real_s(     -18.832407471019791955768596 ) * liVS[52] +
                   li_Real_s(       9.897714159225785124363028 ) * liVS[53] +
                   li_Real_s(      -2.924729950909920717094792 ) * liVS[54] +
                   li_Real_s(       0.374508103554489935049787 ) * liVS[55] +
                   li_Real_s(       0.125308642105039780290099 ) * liVS[56] +
                   li_Real_s(      -0.745138890137234000121680 ) * liVS[57] +
                   li_Real_s(       1.906481486019465876324830 ) * liVS[58] +
                   li_Real_s(      -2.727700625765336894801294 ) * liVS[59] +
                   li_Real_s(       2.362962971940248735336354 ) * liVS[60] +
                   li_Real_s(      -1.241898153585778041829712 ) * liVS[61] +
                   li_Real_s(       0.366975310400102117114329 ) * liVS[62] +
                   li_Real_s(      -0.046990740976568190490070 ) * liVS[63];
        liVC[36] = li_Real_s(       0.151234567739336966951669 ) * liVS[ 0] +
                   li_Real_s(      -0.899305554638786475152301 ) * liVS[ 1] +
                   li_Real_s(       2.300925923866721234389843 ) * liVS[ 2] +
                   li_Real_s(      -3.292052466753489170514513 ) * liVS[ 3] +
                   li_Real_s(       2.851851850316123204720498 ) * liVS[ 4] +
                   li_Real_s(      -1.498842592032890053133087 ) * liVS[ 5] +
                   li_Real_s(       0.442901234453830117132611 ) * liVS[ 6] +
                   li_Real_s(      -0.056712962950898671010691 ) * liVS[ 7] +
                   li_Real_s(      -0.899305555226138864099994 ) * liVS[ 8] +
                   li_Real_s(       5.347656249066099931610552 ) * liVS[ 9] +
                   li_Real_s(     -13.682291667915134070199201 ) * liVS[10] +
                   li_Real_s(      19.575954868946602971391258 ) * liVS[11] +
                   li_Real_s(     -16.958333345429604577248028 ) * liVS[12] +
                   li_Real_s(       8.912760425486638382608362 ) * liVS[13] +
                   li_Real_s(      -2.633680558703323626446036 ) * liVS[14] +
                   li_Real_s(       0.337239583774866957810445 ) * liVS[15] +
                   li_Real_s(       2.300925926548838873486602 ) * liVS[16] +
                   li_Real_s(     -13.682291674694429417513675 ) * liVS[17] +
                   li_Real_s(      35.006944478660145136927895 ) * liVS[18] +
                   li_Real_s(     -50.086226922378372705679794 ) * liVS[19] +
                   li_Real_s(      43.388888968456946315654932 ) * liVS[20] +
                   li_Real_s(     -22.803819494773904352769023 ) * liVS[21] +
                   li_Real_s(       6.738425942697618964416506 ) * liVS[22] +
                   li_Real_s(      -0.862847224516862354448676 ) * liVS[23] +
                   li_Real_s(      -3.292052471791521384147927 ) * liVS[24] +
                   li_Real_s(      19.575954885160673057953318 ) * liVS[25] +
                   li_Real_s(     -50.086226938228541882835998 ) * liVS[26] +
                   li_Real_s(      71.660928016712631460904959 ) * liVS[27] +
                   li_Real_s(     -62.078703876107013570617710 ) * liVS[28] +
                   li_Real_s(      32.626591540562870363828551 ) * liVS[29] +
                   li_Real_s(      -9.641010836907131675843630 ) * liVS[30] +
                   li_Real_s(       1.234519680598041846408819 ) * liVS[31] +
                   li_Real_s(       2.851851855419560877180629 ) * liVS[32] +
                   li_Real_s(     -16.958333363206271826584270 ) * liVS[33] +
                   li_Real_s(      43.388888990609395079900423 ) * liVS[34] +
                   li_Real_s(     -62.078703887487421297919354 ) * liVS[35] +
                   li_Real_s(      53.777777968981183676078217 ) * liVS[36] +
                   li_Real_s(     -28.263889003922493259324256 ) * liVS[37] +
                   li_Real_s(       8.351851889074280776981141 ) * liVS[38] +
                   li_Real_s(      -1.069444449468250013524084 ) * liVS[39] +
                   li_Real_s(      -1.498842595010483336182006 ) * liVS[40] +
                   li_Real_s(       8.912760436248358786315293 ) * liVS[41] +
                   li_Real_s(     -22.803819509431150436284952 ) * liVS[42] +
                   li_Real_s(      32.626591550408846842401545 ) * liVS[43] +
                   li_Real_s(     -28.263889007161907329646056 ) * liVS[44] +
                   li_Real_s(      14.854600764931042533589789 ) * liVS[45] +
                   li_Real_s(      -4.389467615253083820903157 ) * liVS[46] +
                   li_Real_s(       0.562065975268371764705932 ) * liVS[47] +
                   li_Real_s(       0.442901235406226945201524 ) * liVS[48] +
                   li_Real_s(      -2.633680562233186250864492 ) * liVS[49] +
                   li_Real_s(       6.738425947774842938997608 ) * liVS[50] +
                   li_Real_s(      -9.641010840779557611313066 ) * liVS[51] +
                   li_Real_s(       8.351851890847456161282025 ) * liVS[52] +
                   li_Real_s(      -4.389467615687099311116981 ) * liVS[53] +
                   li_Real_s(       1.297067908624352838842242 ) * liVS[54] +
                   li_Real_s(      -0.166087963953032158315182 ) * liVS[55] +
                   li_Real_s(      -0.056712963081281486665830 ) * liVS[56] +
                   li_Real_s(       0.337239584268783421094895 ) * liVS[57] +
                   li_Real_s(      -0.862847225258981609385955 ) * liVS[58] +
                   li_Real_s(       1.234519681215374475868884 ) * liVS[59] +
                   li_Real_s(      -1.069444449799847873805447 ) * liVS[60] +
                   li_Real_s(       0.562065975380224180923960 ) * liVS[61] +
                   li_Real_s(      -0.166087963969859586654820 ) * liVS[62] +
                   li_Real_s(       0.021267361245641325240285 ) * liVS[63];
        liVC[37] = li_Real_s(      -0.024845679016218369383751 ) * liVS[ 0] +
                   li_Real_s(       0.147743055607685302987875 ) * liVS[ 1] +
                   li_Real_s(      -0.378009259504649364203033 ) * liVS[ 2] +
                   li_Real_s(       0.540837191901703207008723 ) * liVS[ 3] +
                   li_Real_s(      -0.468518519161663427041731 ) * liVS[ 4] +
                   li_Real_s(       0.246238426344896232933479 ) * liVS[ 5] +
                   li_Real_s(      -0.072762345820884199998346 ) * liVS[ 6] +
                   li_Real_s(       0.009317129649099920030153 ) * liVS[ 7] +
                   li_Real_s(       0.159336419921106475783290 ) * liVS[ 8] +
                   li_Real_s(      -0.947482640230378647672183 ) * liVS[ 9] +
                   li_Real_s(       2.424189819346774132924338 ) * liVS[10] +
                   li_Real_s(      -3.468412431074543889053530 ) * liVS[11] +
                   li_Real_s(       3.004629638275768677146971 ) * liVS[12] +
                   li_Real_s(      -1.579137736731097918507771 ) * liVS[13] +
                   li_Real_s(       0.466628088132111873420627 ) * liVS[14] +
                   li_Real_s(      -0.059751157639741592220162 ) * liVS[15] +
                   li_Real_s(      -0.437500000772415909011670 ) * liVS[16] +
                   li_Real_s(       2.601562505850772311077890 ) * liVS[17] +
                   li_Real_s(      -6.656250018766766807232216 ) * liVS[18] +
                   li_Real_s(       9.523437532706454433650833 ) * liVS[19] +
                   li_Real_s(      -8.250000033309250468960272 ) * liVS[20] +
                   li_Real_s(       4.335937519813430007786792 ) * liVS[21] +
                   li_Real_s(      -1.281250006384421435967624 ) * liVS[22] +
                   li_Real_s(       0.164062500862198312745477 ) * liVS[23] +
                   li_Real_s(       0.667052470663859864430378 ) * liVS[24] +
                   li_Real_s(      -3.966579872463523059877843 ) * liVS[25] +
                   li_Real_s(      10.148726887523496387188970 ) * liVS[26] +
                   li_Real_s(     -14.520302916014742322659004 ) * liVS[27] +
                   li_Real_s(      12.578703765058826036238315 ) * liVS[28] +
                   li_Real_s(      -6.610966471323675541782450 ) * liVS[29] +
                   li_Real_s(       1.953510814039229437355516 ) * liVS[30] +
                   li_Real_s(      -0.250144677483457589239890 ) * liVS[31] +
                   li_Real_s(      -0.610339507806767045394736 ) * liVS[32] +
                   li_Real_s(       3.629340289796019636270330 ) * liVS[33] +
                   li_Real_s(      -9.285879666971636581251914 ) * liVS[34] +
                   li_Real_s(      13.285783242309429752481265 ) * liVS[35] +
                   li_Real_s(     -11.509259322331592656496468 ) * liVS[36] +
                   li_Real_s(       6.048900499884314996279500 ) * liVS[37] +
                   li_Real_s(      -1.787422851274402146870557 ) * liVS[38] +
                   li_Real_s(       0.228877316394634711116396 ) * liVS[39] +
                   li_Real_s(       0.335416667659891976072117 ) * liVS[40] +
                   li_Real_s(      -1.994531257263061352347222 ) * liVS[41] +
                   li_Real_s(       5.103125022408530497841639 ) * liVS[42] +
                   li_Real_s(      -7.301302121068554917826532 ) * liVS[43] +
                   li_Real_s(       6.325000037394267415891136 ) * liVS[44] +
                   li_Real_s(      -3.324218771792746274229557 ) * liVS[45] +
                   li_Real_s(       0.982291673588533420513613 ) * liVS[46] +
                   li_Real_s(      -0.125781250926860793670770 ) * liVS[47] +
                   li_Real_s(      -0.102623457114470628681602 ) * liVS[48] +
                   li_Real_s(       0.610243057919886755513517 ) * liVS[49] +
                   li_Real_s(      -1.561342599854685708749003 ) * liVS[50] +
                   li_Real_s(       2.233892759088862867145053 ) * liVS[51] +
                   li_Real_s(      -1.935185197203603379989545 ) * liVS[52] +
                   li_Real_s(       1.017071766240105201717370 ) * liVS[53] +
                   li_Real_s(      -0.300540125667767554773491 ) * liVS[54] +
                   li_Real_s(       0.038483796591672003728490 ) * liVS[55] +
                   li_Real_s(       0.013503086464071500927275 ) * liVS[56] +
                   li_Real_s(      -0.080295139211537955303832 ) * liVS[57] +
                   li_Real_s(       0.205439815803164615992671 ) * liVS[58] +
                   li_Real_s(      -0.293933257825042204558486 ) * liVS[59] +
                   li_Real_s(       0.254629631256096278235646 ) * liVS[60] +
                   li_Real_s(      -0.133825232423831541606063 ) * liVS[61] +
                   li_Real_s(       0.039544753384181341449022 ) * liVS[62] +
                   li_Real_s(      -0.005063657447100133879303 ) * liVS[63];
        liVC[38] = li_Real_s(       0.002160493831729004554632 ) * liVS[ 0] +
                   li_Real_s(      -0.012847222254476764469189 ) * liVS[ 1] +
                   li_Real_s(       0.032870370470016896380372 ) * liVS[ 2] +
                   li_Real_s(      -0.047029321158175485528830 ) * liVS[ 3] +
                   li_Real_s(       0.040740740912966366682468 ) * liVS[ 4] +
                   li_Real_s(      -0.021412037139119211026639 ) * liVS[ 5] +
                   li_Real_s(       0.006327160526666923689998 ) * liVS[ 6] +
                   li_Real_s(      -0.000810185189609916034392 ) * liVS[ 7] +
                   li_Real_s(      -0.014583333377391538565604 ) * liVS[ 8] +
                   li_Real_s(       0.086718750310998224395398 ) * liVS[ 9] +
                   li_Real_s(      -0.221875000946187950212618 ) * liVS[10] +
                   li_Real_s(       0.317447918254831762396861 ) * liVS[11] +
                   li_Real_s(      -0.275000001576674646397436 ) * liVS[12] +
                   li_Real_s(       0.144531250922804499481344 ) * liVS[13] +
                   li_Real_s(      -0.042708333628033973150551 ) * liVS[14] +
                   li_Real_s(       0.005468750039653635930392 ) * liVS[15] +
                   li_Real_s(       0.042129629782887367817068 ) * liVS[16] +
                   li_Real_s(      -0.250520834412952864145296 ) * liVS[17] +
                   li_Real_s(       0.640972225480433532851521 ) * liVS[18] +
                   li_Real_s(      -0.917071764676351763867501 ) * liVS[19] +
                   li_Real_s(       0.794444449776718464306668 ) * liVS[20] +
                   li_Real_s(      -0.417534725322769917355004 ) * liVS[21] +
                   li_Real_s(       0.123379630615600754950378 ) * liVS[22] +
                   li_Real_s(      -0.015798611243565519046683 ) * liVS[23] +
                   li_Real_s(      -0.067515432369249728239424 ) * liVS[24] +
                   li_Real_s(       0.401475696346249188994904 ) * liVS[25] +
                   li_Real_s(      -1.027199079784627766542826 ) * liVS[26] +
                   li_Real_s(       1.469666290304140643385722 ) * liVS[27] +
                   li_Real_s(      -1.273148157392013368749417 ) * liVS[28] +
                   li_Real_s(       0.669126162759850373618065 ) * liVS[29] +
                   li_Real_s(      -0.197723767128948402227451 ) * liVS[30] +
                   li_Real_s(       0.025318287264599677321986 ) * liVS[31] +
                   li_Real_s(       0.064814815087421262518319 ) * liVS[32] +
                   li_Real_s(      -0.385416668580115584497037 ) * liVS[33] +
                   li_Real_s(       0.986111116835310497208411 ) * liVS[34] +
                   li_Real_s(      -1.410879639053096434864187 ) * liVS[35] +
                   li_Real_s(       1.222222231414889836997872 ) * liVS[36] +
                   li_Real_s(      -0.642361116416698130215934 ) * liVS[37] +
                   li_Real_s(       0.189814816492407639758611 ) * liVS[38] +
                   li_Real_s(      -0.024305555780118948128177 ) * liVS[39] +
                   li_Real_s(      -0.037268518678583606451582 ) * liVS[40] +
                   li_Real_s(       0.221614584455389729278352 ) * liVS[41] +
                   li_Real_s(      -0.567013892236121863632548 ) * liVS[42] +
                   li_Real_s(       0.811255792530002017848290 ) * liVS[43] +
                   li_Real_s(      -0.702777783120207888423181 ) * liVS[44] +
                   li_Real_s(       0.369357641964068339479610 ) * liVS[45] +
                   li_Real_s(      -0.109143519488601181688381 ) * liVS[46] +
                   li_Real_s(       0.013975694574054344301861 ) * liVS[47] +
                   li_Real_s(       0.011882716100563628369713 ) * liVS[48] +
                   li_Real_s(      -0.070659722580763750787725 ) * liVS[49] +
                   li_Real_s(       0.180787038104494479640039 ) * liVS[50] +
                   li_Real_s(      -0.258661267179747467181983 ) * liVS[51] +
                   li_Real_s(       0.224074075769910896482884 ) * liVS[52] +
                   li_Real_s(      -0.117766204677731634031801 ) * liVS[53] +
                   li_Real_s(       0.034799383022687591093813 ) * liVS[54] +
                   li_Real_s(      -0.004456018559413743584940 ) * liVS[55] +
                   li_Real_s(      -0.001620370377282021046028 ) * liVS[56] +
                   li_Real_s(       0.009635416715073025911709 ) * liVS[57] +
                   li_Real_s(      -0.024652777921705393282537 ) * liVS[58] +
                   li_Real_s(       0.035271990975983602556454 ) * liVS[59] +
                   li_Real_s(      -0.030555555783424614979538 ) * liVS[60] +
                   li_Real_s(       0.016059027908427683795090 ) * liVS[61] +
                   li_Real_s(      -0.004745370411428403989440 ) * liVS[62] +
                   li_Real_s(       0.000607638894355999603647 ) * liVS[63];
        liVC[39] = li_Real_s(      -0.000077160494183173278238 ) * liVS[ 0] +
                   li_Real_s(       0.000458829367481539129869 ) * liVS[ 1] +
                   li_Real_s(      -0.001173941805999671045280 ) * liVS[ 2] +
                   li_Real_s(       0.001679618618251631223526 ) * liVS[ 3] +
                   li_Real_s(      -0.001455026466292287656845 ) * liVS[ 4] +
                   li_Real_s(       0.000764715614983115660186 ) * liVS[ 5] +
                   li_Real_s(      -0.000225970019703328629107 ) * liVS[ 6] +
                   li_Real_s(       0.000028935185462132528844 ) * liVS[ 7] +
                   li_Real_s(       0.000540123459681730172521 ) * liVS[ 8] +
                   li_Real_s(      -0.003211805575274000613267 ) * liVS[ 9] +
                   li_Real_s(       0.008217592650610375595854 ) * liVS[10] +
                   li_Real_s(      -0.011757330341569435794291 ) * liVS[11] +
                   li_Real_s(       0.010185185277124144190286 ) * liVS[12] +
                   li_Real_s(      -0.005353009312255060936359 ) * liVS[13] +
                   li_Real_s(       0.001581790140226531360401 ) * liVS[14] +
                   li_Real_s(      -0.000202546298544284408827 ) * liVS[15] +
                   li_Real_s(      -0.001620370379655414194708 ) * liVS[16] +
                   li_Real_s(       0.009635416730225060052972 ) * liVS[17] +
                   li_Real_s(      -0.024652777964664186055588 ) * liVS[18] +
                   li_Real_s(       0.035271991044837233686415 ) * liVS[19] +
                   li_Real_s(      -0.030555555850011593965743 ) * liVS[20] +
                   li_Real_s(       0.016059027947050344420354 ) * liVS[21] +
                   li_Real_s(      -0.004745370423832268383380 ) * liVS[22] +
                   li_Real_s(       0.000607638896050817500782 ) * liVS[23] +
                   li_Real_s(       0.002700617299662214088052 ) * liVS[24] +
                   li_Real_s(      -0.016059027885469787794159 ) * liVS[25] +
                   li_Real_s(       0.041087963279219673340403 ) * liVS[26] +
                   li_Real_s(      -0.058786651747984794424973 ) * liVS[27] +
                   li_Real_s(       0.050925926421832939938739 ) * liVS[28] +
                   li_Real_s(      -0.026765046580734842829319 ) * liVS[29] +
                   li_Real_s(       0.007908950706955058390646 ) * liVS[30] +
                   li_Real_s(      -0.001012731493480466780921 ) * liVS[31] +
                   li_Real_s(      -0.002700617299402685578258 ) * liVS[32] +
                   li_Real_s(       0.016059027883709574291426 ) * liVS[33] +
                   li_Real_s(      -0.041087963273576111333796 ) * liVS[34] +
                   li_Real_s(       0.058786651737738893075402 ) * liVS[35] +
                   li_Real_s(      -0.050925926410856761639145 ) * liVS[36] +
                   li_Real_s(       0.026765046573862118217679 ) * liVS[37] +
                   li_Real_s(      -0.007908950704627206937380 ) * liVS[38] +
                   li_Real_s(       0.001012731493152178169348 ) * liVS[39] +
                   li_Real_s(       0.001620370379306552630072 ) * liVS[40] +
                   li_Real_s(      -0.009635416727917747725662 ) * liVS[41] +
                   li_Real_s(       0.024652777957114679896478 ) * liVS[42] +
                   li_Real_s(      -0.035271991030696614455309 ) * liVS[43] +
                   li_Real_s(       0.030555555834448609037191 ) * liVS[44] +
                   li_Real_s(      -0.016059027937107572120423 ) * liVS[45] +
                   li_Real_s(       0.004745370420417888493447 ) * liVS[46] +
                   li_Real_s(      -0.000607638895565795755793 ) * liVS[47] +
                   li_Real_s(      -0.000540123459620546475524 ) * liVS[48] +
                   li_Real_s(       0.003211805574954569933444 ) * liVS[49] +
                   li_Real_s(      -0.008217592649325602960353 ) * liVS[50] +
                   li_Real_s(       0.011757330338499231980887 ) * liVS[51] +
                   li_Real_s(      -0.010185185273132495165083 ) * liVS[52] +
                   li_Real_s(       0.005353009309420514102995 ) * liVS[53] +
                   li_Real_s(      -0.001581790139187449345526 ) * liVS[54] +
                   li_Real_s(       0.000202546298391781398607 ) * liVS[55] +
                   li_Real_s(       0.000077160494207369201280 ) * liVS[56] +
                   li_Real_s(      -0.000458829367685039974517 ) * liVS[57] +
                   li_Real_s(       0.001173941806555717573546 ) * liVS[58] +
                   li_Real_s(      -0.001679618618978639954520 ) * liVS[59] +
                   li_Real_s(       0.001455026466799898970939 ) * liVS[60] +
                   li_Real_s(      -0.000764715615171341722527 ) * liVS[61] +
                   li_Real_s(       0.000225970019736562461460 ) * liVS[62] +
                   li_Real_s(      -0.000028935185464516038900 ) * liVS[63];
        liVC[40] = li_Real_s(      -0.063888888935330567786508 ) * liVS[ 0] +
                   li_Real_s(       0.409722222529569080062828 ) * liVS[ 1] +
                   li_Real_s(      -1.125000000876775096969595 ) * liVS[ 2] +
                   li_Real_s(       1.715277779167660110459792 ) * liVS[ 3] +
                   li_Real_s(      -1.569444445763851225805752 ) * liVS[ 4] +
                   li_Real_s(       0.862500000750080375588880 ) * liVS[ 5] +
                   li_Real_s(      -0.263888889125574288563314 ) * liVS[ 6] +
                   li_Real_s(       0.034722222254222286086378 ) * liVS[ 7] +
                   li_Real_s(       0.000000000288088293053672 ) * liVS[ 8] +
                   li_Real_s(      -0.000000001913550583065704 ) * liVS[ 9] +
                   li_Real_s(       0.000000005469710313620106 ) * liVS[10] +
                   li_Real_s(      -0.000000008677916515220360 ) * liVS[11] +
                   li_Real_s(       0.000000008238646470365740 ) * liVS[12] +
                   li_Real_s(      -0.000000004682127537088978 ) * liVS[13] +
                   li_Real_s(       0.000000001476688300047327 ) * liVS[14] +
                   li_Real_s(      -0.000000000199538741711805 ) * liVS[15] +
                   li_Real_s(      -0.000000000766885163365755 ) * liVS[16] +
                   li_Real_s(       0.000000005112637011911346 ) * liVS[17] +
                   li_Real_s(      -0.000000014643506772074770 ) * liVS[18] +
                   li_Real_s(       0.000000023251727048854645 ) * liVS[19] +
                   li_Real_s(      -0.000000022075831577627204 ) * liVS[20] +
                   li_Real_s(       0.000000012541026048181057 ) * liVS[21] +
                   li_Real_s(      -0.000000003952977963192395 ) * liVS[22] +
                   li_Real_s(       0.000000000533811367313079 ) * liVS[23] +
                   li_Real_s(       0.000000001131527030851618 ) * liVS[24] +
                   li_Real_s(      -0.000000007573192490198220 ) * liVS[25] +
                   li_Real_s(       0.000000021739994065963193 ) * liVS[26] +
                   li_Real_s(      -0.000000034555388841784584 ) * liVS[27] +
                   li_Real_s(       0.000000032814324757485434 ) * liVS[28] +
                   li_Real_s(      -0.000000018635979770802171 ) * liVS[29] +
                   li_Real_s(       0.000000005871074944478038 ) * liVS[30] +
                   li_Real_s(      -0.000000000792359695993281 ) * liVS[31] +
                   li_Real_s(      -0.000000000998729753096154 ) * liVS[32] +
                   li_Real_s(       0.000000006712659913411519 ) * liVS[33] +
                   li_Real_s(      -0.000000019319505715529515 ) * liVS[34] +
                   li_Real_s(       0.000000030748519060840803 ) * liVS[35] +
                   li_Real_s(      -0.000000029211981166957491 ) * liVS[36] +
                   li_Real_s(       0.000000016588395410646890 ) * liVS[37] +
                   li_Real_s(      -0.000000005224052656987978 ) * liVS[38] +
                   li_Real_s(       0.000000000704694907671930 ) * liVS[39] +
                   li_Real_s(       0.000000000527895523638871 ) * liVS[40] +
                   li_Real_s(      -0.000000003564088727231579 ) * liVS[41] +
                   li_Real_s(       0.000000010287369052600213 ) * liVS[42] +
                   li_Real_s(      -0.000000016399446942509929 ) * liVS[43] +
                   li_Real_s(       0.000000015590608746123144 ) * liVS[44] +
                   li_Real_s(      -0.000000008854236188456930 ) * liVS[45] +
                   li_Real_s(       0.000000002787827630829214 ) * liVS[46] +
                   li_Real_s(      -0.000000000375929094993012 ) * liVS[47] +
                   li_Real_s(      -0.000000000155083102649231 ) * liVS[48] +
                   li_Real_s(       0.000000001051902187478170 ) * liVS[49] +
                   li_Real_s(      -0.000000003045565718712260 ) * liVS[50] +
                   li_Real_s(       0.000000004863813130145583 ) * liVS[51] +
                   li_Real_s(      -0.000000004627964749867781 ) * liVS[52] +
                   li_Real_s(       0.000000002629032149397634 ) * liVS[53] +
                   li_Real_s(      -0.000000000827723726708332 ) * liVS[54] +
                   li_Real_s(       0.000000000111589830916217 ) * liVS[55] +
                   li_Real_s(       0.000000000019600231042101 ) * liVS[56] +
                   li_Real_s(      -0.000000000133547638766100 ) * liVS[57] +
                   li_Real_s(       0.000000000387850267628128 ) * liVS[58] +
                   li_Real_s(      -0.000000000620554616034307 ) * liVS[59] +
                   li_Real_s(       0.000000000591027780100234 ) * liVS[60] +
                   li_Real_s(      -0.000000000335869460133700 ) * liVS[61] +
                   li_Real_s(       0.000000000105747975019474 ) * liVS[62] +
                   li_Real_s(      -0.000000000014254538855856 ) * liVS[63];
        liVC[41] = li_Real_s(       0.165654761859979515747909 ) * liVS[ 0] +
                   li_Real_s(      -1.062351190232186937123515 ) * liVS[ 1] +
                   li_Real_s(       2.916964285209871388104830 ) * liVS[ 2] +
                   li_Real_s(      -4.447470237601013565154062 ) * liVS[ 3] +
                   li_Real_s(       4.069345237882513011129504 ) * liVS[ 4] +
                   li_Real_s(      -2.236339285710883473257127 ) * liVS[ 5] +
                   li_Real_s(       0.684226190502207787247357 ) * liVS[ 6] +
                   li_Real_s(      -0.090029761910498606880537 ) * liVS[ 7] +
                   li_Real_s(      -0.447222222126029578248563 ) * liVS[ 8] +
                   li_Real_s(       2.868055555318298033284918 ) * liVS[ 9] +
                   li_Real_s(      -7.875000000615459683217523 ) * liVS[10] +
                   li_Real_s(      12.006944447372648099303660 ) * liVS[11] +
                   li_Real_s(     -10.986111115476765220932975 ) * liVS[12] +
                   li_Real_s(       6.037500003185494179547277 ) * liVS[13] +
                   li_Real_s(      -1.847222223382288142090601 ) * liVS[14] +
                   li_Real_s(       0.243055555724103200532227 ) * liVS[15] +
                   li_Real_s(       0.670833333492392824837225 ) * liVS[16] +
                   li_Real_s(      -4.302083335614064196761319 ) * liVS[17] +
                   li_Real_s(      11.812500010163873298552062 ) * liVS[18] +
                   li_Real_s(     -18.010416688083239478146425 ) * liVS[19] +
                   li_Real_s(      16.479166691184840232153874 ) * liVS[20] +
                   li_Real_s(      -9.056250015738914527219094 ) * liVS[21] +
                   li_Real_s(       2.770833338676506052422610 ) * liVS[22] +
                   li_Real_s(      -0.364583334081395094017353 ) * liVS[23] +
                   li_Real_s(      -0.745370371090634975530520 ) * liVS[24] +
                   li_Real_s(       4.780092599395806551854093 ) * liVS[25] +
                   li_Real_s(     -13.125000025116065671682009 ) * liVS[26] +
                   li_Real_s(      20.011574121897766076472180 ) * liVS[27] +
                   li_Real_s(     -18.310185236751241433239556 ) * liVS[28] +
                   li_Real_s(      10.062500031888678364566658 ) * liVS[29] +
                   li_Real_s(      -3.078703714276793235171681 ) * liVS[30] +
                   li_Real_s(       0.405092594052494092693451 ) * liVS[31] +
                   li_Real_s(       0.559027778757007709486970 ) * liVS[32] +
                   li_Real_s(      -3.585069452913682397365847 ) * liVS[33] +
                   li_Real_s(       9.843750029480716534635576 ) * liVS[34] +
                   li_Real_s(     -15.008680609549029938420972 ) * liVS[35] +
                   li_Real_s(      13.732638945610828429266803 ) * liVS[36] +
                   li_Real_s(      -7.546875034461829834242508 ) * liVS[37] +
                   li_Real_s(       2.309027789067423341862195 ) * liVS[38] +
                   li_Real_s(      -0.303819445991431180686959 ) * liVS[39] +
                   li_Real_s(      -0.268333334001431111914826 ) * liVS[40] +
                   li_Real_s(       1.720833338906078147090284 ) * liVS[41] +
                   li_Real_s(      -4.725000018864832895815198 ) * liVS[42] +
                   li_Real_s(       7.204166700524773148117674 ) * liVS[43] +
                   li_Real_s(      -6.591666701723187138384219 ) * liVS[44] +
                   li_Real_s(       3.622500021077569964944587 ) * liVS[45] +
                   li_Real_s(      -1.108333340186579807351563 ) * liVS[46] +
                   li_Real_s(       0.145833334267604058931411 ) * liVS[47] +
                   li_Real_s(       0.074537037270374195685463 ) * liVS[48] +
                   li_Real_s(      -0.478009261170884935943093 ) * liVS[49] +
                   li_Real_s(       1.312500006371976724039996 ) * liVS[50] +
                   li_Real_s(      -2.001157418707064650220673 ) * liVS[51] +
                   li_Real_s(       1.831018530112402231679880 ) * liVS[52] +
                   li_Real_s(      -1.006250006923100404065963 ) * liVS[53] +
                   li_Real_s(       0.307870372609668763175250 ) * liVS[54] +
                   li_Real_s(      -0.040509259563370036971719 ) * liVS[55] +
                   li_Real_s(      -0.009126984160230833253991 ) * liVS[56] +
                   li_Real_s(       0.058531746301701520485139 ) * liVS[57] +
                   li_Real_s(      -0.160714286606217005015651 ) * liVS[58] +
                   li_Real_s(       0.245039684109766398023567 ) * liVS[59] +
                   li_Real_s(      -0.224206350807890864018646 ) * liVS[60] +
                   li_Real_s(       0.123214286666148808979671 ) * liVS[61] +
                   li_Real_s(      -0.037698413005132769271199 ) * liVS[62] +
                   li_Real_s(       0.004960317501851108090705 ) * liVS[63];
        liVC[42] = li_Real_s(      -0.166466049134861293623544 ) * liVS[ 0] +
                   li_Real_s(       1.067554010838920319770295 ) * liVS[ 1] +
                   li_Real_s(      -2.931249996172915039949203 ) * liVS[ 2] +
                   li_Real_s(       4.469251537883017988406209 ) * liVS[ 3] +
                   li_Real_s(      -4.089274686917573120581437 ) * liVS[ 4] +
                   li_Real_s(       2.247291664422611745521863 ) * liVS[ 5] +
                   li_Real_s(      -0.687577159849693808624238 ) * liVS[ 6] +
                   li_Real_s(       0.090470678930436587705799 ) * liVS[ 7] +
                   li_Real_s(       0.712361109984181339882525 ) * liVS[ 8] +
                   li_Real_s(      -4.568402771308955223616977 ) * liVS[ 9] +
                   li_Real_s(      12.543749984692372123618043 ) * liVS[10] +
                   li_Real_s(     -19.125347202682888791969162 ) * liVS[11] +
                   li_Real_s(      17.499305540836779471192131 ) * liVS[12] +
                   li_Real_s(      -9.616874993326650411518131 ) * liVS[13] +
                   li_Real_s(       2.942361109376305350338043 ) * liVS[14] +
                   li_Real_s(      -0.387152777571142081569633 ) * liVS[15] +
                   li_Real_s(      -1.403958331249096858073244 ) * liVS[16] +
                   li_Real_s(       9.003645822514922514301361 ) * liVS[17] +
                   li_Real_s(     -24.721874977973985920698397 ) * liVS[18] +
                   li_Real_s(      37.693229144591363422023278 ) * liVS[19] +
                   li_Real_s(     -34.488541655880048608651123 ) * liVS[20] +
                   li_Real_s(      18.953437498250440995661847 ) * liVS[21] +
                   li_Real_s(      -5.798958333699395062410531 ) * liVS[22] +
                   li_Real_s(       0.763020833445796853311549 ) * liVS[23] +
                   li_Real_s(       1.684182096733028544122135 ) * liVS[24] +
                   li_Real_s(     -10.800733016080940274150635 ) * liVS[25] +
                   li_Real_s(      29.656249989108381726055086 ) * liVS[26] +
                   li_Real_s(     -45.216628088504904781075311 ) * liVS[27] +
                   li_Real_s(      41.372299399724973056891031 ) * liVS[28] +
                   li_Real_s(     -22.736458349323509509076757 ) * liVS[29] +
                   li_Real_s(       6.956404327298342238350415 ) * liVS[30] +
                   li_Real_s(      -0.915316358955339470782064 ) * liVS[31] +
                   li_Real_s(      -1.309722221143672982179851 ) * liVS[32] +
                   li_Real_s(       8.399305553063527440826874 ) * liVS[33] +
                   li_Real_s(     -23.062500005508393741138207 ) * liVS[34] +
                   li_Real_s(      35.163194470644910438750230 ) * liVS[35] +
                   li_Real_s(     -32.173611149315178181495867 ) * liVS[36] +
                   li_Real_s(      17.681250027193655682822282 ) * liVS[37] +
                   li_Real_s(      -5.409722231846706641533729 ) * liVS[38] +
                   li_Real_s(       0.711805556911851766699328 ) * liVS[39] +
                   li_Real_s(       0.642083333065748718126997 ) * liVS[40] +
                   li_Real_s(      -4.117708334210934850716512 ) * liVS[41] +
                   li_Real_s(      11.306250009467795791806566 ) * liVS[42] +
                   li_Real_s(     -17.238541691132109434647646 ) * liVS[43] +
                   li_Real_s(      15.772916696991856611020921 ) * liVS[44] +
                   li_Real_s(      -8.668125020086581145051241 ) * liVS[45] +
                   li_Real_s(       2.652083340186157656148680 ) * liVS[46] +
                   li_Real_s(      -0.348958334281927851083793 ) * liVS[47] +
                   li_Real_s(      -0.180841049375469964388685 ) * liVS[48] +
                   li_Real_s(       1.159741513131725465513000 ) * liVS[49] +
                   li_Real_s(      -3.184375004376811091333366 ) * liVS[50] +
                   li_Real_s(       4.855189052997150866985976 ) * liVS[51] +
                   li_Real_s(      -4.442399702755409407473053 ) * liVS[52] +
                   li_Real_s(       2.441354173970770702339905 ) * liVS[53] +
                   li_Real_s(      -0.746952162939390262863526 ) * liVS[54] +
                   li_Real_s(       0.098283179347433247130539 ) * liVS[55] +
                   li_Real_s(       0.022361111117024989880520 ) * liVS[56] +
                   li_Real_s(      -0.143402777928990698974587 ) * liVS[57] +
                   li_Real_s(       0.393750000711876602110806 ) * liVS[58] +
                   li_Real_s(      -0.600347223719616351900186 ) * liVS[59] +
                   li_Real_s(       0.549305557245904019225691 ) * liVS[60] +
                   li_Real_s(      -0.301875001063852677063437 ) * liVS[61] +
                   li_Real_s(       0.092361111463343803507087 ) * liVS[62] +
                   li_Real_s(      -0.012152777825697569369368 ) * liVS[63];
        liVC[43] = li_Real_s(       0.085806326984503300536744 ) * liVS[ 0] +
                   li_Real_s(      -0.550279705736711477470635 ) * liVS[ 1] +
                   li_Real_s(       1.510937497373938676048510 ) * liVS[ 2] +
                   li_Real_s(      -2.303713345190427475017714 ) * liVS[ 3] +
                   li_Real_s(       2.107851077338256828852536 ) * liVS[ 4] +
                   li_Real_s(      -1.158385415230238635331261 ) * liVS[ 5] +
                   li_Real_s(       0.354417437866284434910824 ) * liVS[ 6] +
                   li_Real_s(      -0.046633873405570902548334 ) * liVS[ 7] +
                   li_Real_s(      -0.452901233815477866073707 ) * liVS[ 8] +
                   li_Real_s(       2.904475304492611265061441 ) * liVS[ 9] +
                   li_Real_s(      -7.974999990717051900901424 ) * liVS[10] +
                   li_Real_s(      12.159413569318083148118603 ) * liVS[11] +
                   li_Real_s(     -11.125617276613112949235074 ) * liVS[12] +
                   li_Real_s(       6.114166663815926661129652 ) * liVS[13] +
                   li_Real_s(      -1.870679011724656781723297 ) * liVS[14] +
                   li_Real_s(       0.246141975243677535445386 ) * liVS[15] +
                   li_Real_s(       1.045914350600412490166491 ) * liVS[16] +
                   li_Real_s(      -6.707494207185195733700311 ) * liVS[17] +
                   li_Real_s(      18.417187490731876664540323 ) * liVS[18] +
                   li_Real_s(     -28.080526616024592101439339 ) * liVS[19] +
                   li_Real_s(      25.693113429908962075387535 ) * liVS[20] +
                   li_Real_s(     -14.119843755695015374840295 ) * liVS[21] +
                   li_Real_s(       4.320081021004123833506583 ) * liVS[22] +
                   li_Real_s(      -0.568431713340576294513085 ) * liVS[23] +
                   li_Real_s(      -1.380709875553868926090217 ) * liVS[24] +
                   li_Real_s(       8.854552466662656939888620 ) * liVS[25] +
                   li_Real_s(     -24.312500004134914632913933 ) * liVS[26] +
                   li_Real_s(      37.069058664166121275229671 ) * liVS[27] +
                   li_Real_s(     -33.917438304686761796347128 ) * liVS[28] +
                   li_Real_s(      18.639583357094632987127625 ) * liVS[29] +
                   li_Real_s(      -5.702932107206748923999839 ) * liVS[30] +
                   li_Real_s(       0.750385803658884631417436 ) * liVS[31] +
                   li_Real_s(       1.129147376262221058595969 ) * liVS[32] +
                   li_Real_s(      -7.241271221180503481207325 ) * liVS[33] +
                   li_Real_s(      19.882812516163095750698631 ) * liVS[34] +
                   li_Real_s(     -30.315152431721202219705447 ) * liVS[35] +
                   li_Real_s(      27.737750819998538531763188 ) * liVS[36] +
                   li_Real_s(     -15.243489615129993453024326 ) * liVS[37] +
                   li_Real_s(       4.663869609568883412009654 ) * liVS[38] +
                   li_Real_s(      -0.613667053961039599130345 ) * liVS[39] +
                   li_Real_s(      -0.570740740844994931535439 ) * liVS[40] +
                   li_Real_s(       3.660185188096113595435099 ) * liVS[41] +
                   li_Real_s(     -10.050000014001307491184889 ) * liVS[42] +
                   li_Real_s(      15.323148177984780460292313 ) * liVS[43] +
                   li_Real_s(     -14.020370404366788719130454 ) * liVS[44] +
                   li_Real_s(       7.705000021547382615949573 ) * liVS[45] +
                   li_Real_s(      -2.357407414580461768593977 ) * liVS[46] +
                   li_Real_s(       0.310185186165275794678564 ) * liVS[47] +
                   li_Real_s(       0.164070216140036251317724 ) * liVS[48] +
                   li_Real_s(      -1.052189430310472850749193 ) * liVS[49] +
                   li_Real_s(       2.889062505447714812589766 ) * liVS[50] +
                   li_Real_s(      -4.404928637479002873078571 ) * liVS[51] +
                   li_Real_s(       4.030420536752899351995438 ) * liVS[52] +
                   li_Real_s(      -2.214947924164150094838988 ) * liVS[53] +
                   li_Real_s(       0.677681329626594219917024 ) * liVS[54] +
                   li_Real_s(      -0.089168596013619705331621 ) * liVS[55] +
                   li_Real_s(      -0.020586419770292962994063 ) * liVS[56] +
                   li_Real_s(       0.132021605146006359987609 ) * liVS[57] +
                   li_Real_s(      -0.362500000821680323781493 ) * liVS[58] +
                   li_Real_s(       0.552700618884001571018416 ) * liVS[59] +
                   li_Real_s(      -0.505709878276272561947735 ) * liVS[60] +
                   li_Real_s(       0.277916667731450905964863 ) * liVS[61] +
                   li_Real_s(      -0.085030864545030393486513 ) * liVS[62] +
                   li_Real_s(       0.011188271651828785024918 ) * liVS[63];
        liVC[44] = li_Real_s(      -0.024845678974488194512560 ) * liVS[ 0] +
                   li_Real_s(       0.159336419543855356550921 ) * liVS[ 1] +
                   li_Real_s(      -0.437499999540319262081312 ) * liVS[ 2] +
                   li_Real_s(       0.667052468618718918946797 ) * liVS[ 3] +
                   li_Real_s(      -0.610339505855495012554002 ) * liVS[ 4] +
                   li_Real_s(       0.335416666562532661544083 ) * liVS[ 5] +
                   li_Real_s(      -0.102623456773322074830901 ) * liVS[ 6] +
                   li_Real_s(       0.013503086418520160449930 ) * liVS[ 7] +
                   li_Real_s(       0.147743055458153360604001 ) * liVS[ 8] +
                   li_Real_s(      -0.947482638560153800355579 ) * liVS[ 9] +
                   li_Real_s(       2.601562500004292122213201 ) * liVS[10] +
                   li_Real_s(      -3.966579862499120423535715 ) * liVS[11] +
                   li_Real_s(       3.629340280192504764045225 ) * liVS[12] +
                   li_Real_s(      -1.994531251839496865940760 ) * liVS[13] +
                   li_Real_s(       0.610243056228978675292751 ) * liVS[14] +
                   li_Real_s(      -0.080295138985157166189310 ) * liVS[15] +
                   li_Real_s(      -0.378009259306328004868192 ) * liVS[16] +
                   li_Real_s(       2.424189816092368499056420 ) * liVS[17] +
                   li_Real_s(      -6.656250006334894386839096 ) * liVS[18] +
                   li_Real_s(      10.148726865697767607343849 ) * liVS[19] +
                   li_Real_s(      -9.285879645718477348736997 ) * liVS[20] +
                   li_Real_s(       5.103125010360973767831183 ) * liVS[21] +
                   li_Real_s(      -1.561342596089365653000414 ) * liVS[22] +
                   li_Real_s(       0.205439815297955519213247 ) * liVS[23] +
                   li_Real_s(       0.540837191785026760726396 ) * liVS[24] +
                   li_Real_s(      -3.468412427247352880499420 ) * liVS[25] +
                   li_Real_s(       9.523437516770226096696206 ) * liVS[26] +
                   li_Real_s(     -14.520302887291203219888303 ) * liVS[27] +
                   li_Real_s(      13.285783214095118864861433 ) * liVS[28] +
                   li_Real_s(      -7.301302105026699962309067 ) * liVS[29] +
                   li_Real_s(       2.233892754064483732179269 ) * liVS[30] +
                   li_Real_s(      -0.293933257149592064294552 ) * liVS[31] +
                   li_Real_s(      -0.468518519151102097453077 ) * liVS[32] +
                   li_Real_s(       3.004629635336039772397498 ) * liVS[33] +
                   li_Real_s(      -8.250000020159552249765511 ) * liVS[34] +
                   li_Real_s(      12.578703740872928307226175 ) * liVS[35] +
                   li_Real_s(     -11.509259298414704630886263 ) * liVS[36] +
                   li_Real_s(       6.325000023761788270348916 ) * liVS[37] +
                   li_Real_s(      -1.935185192925262587237967 ) * liVS[38] +
                   li_Real_s(       0.254629630679866103548648 ) * liVS[39] +
                   li_Real_s(       0.246238426373746932540598 ) * liVS[40] +
                   li_Real_s(      -1.579137735310591539317215 ) * liVS[41] +
                   li_Real_s(       4.335937513064070714108311 ) * liVS[42] +
                   li_Real_s(      -6.610966458710437976264984 ) * liVS[43] +
                   li_Real_s(       6.048900487346949184086498 ) * liVS[44] +
                   li_Real_s(      -3.324218764632895783961430 ) * liVS[45] +
                   li_Real_s(       1.017071763989368271552394 ) * liVS[46] +
                   li_Real_s(      -0.133825232120211551345434 ) * liVS[47] +
                   li_Real_s(      -0.072762345838578657009066 ) * liVS[48] +
                   li_Real_s(       0.466628087748009567903296 ) * liVS[49] +
                   li_Real_s(      -1.281250004444503431955127 ) * liVS[50] +
                   li_Real_s(       1.953510810360667093732445 ) * liVS[51] +
                   li_Real_s(      -1.787422847602684328194300 ) * liVS[52] +
                   li_Real_s(       0.982291671489348150458909 ) * liVS[53] +
                   li_Real_s(      -0.300540125007238145826705 ) * liVS[54] +
                   li_Real_s(       0.039544753294980417024362 ) * liVS[55] +
                   li_Real_s(       0.009317129652609779100203 ) * liVS[56] +
                   li_Real_s(      -0.059751157596096948765307 ) * liVS[57] +
                   li_Real_s(       0.164062500624307716634576 ) * liVS[58] +
                   li_Real_s(      -0.250144677024831452172293 ) * liVS[59] +
                   li_Real_s(       0.228877315934798541974260 ) * liVS[60] +
                   li_Real_s(      -0.125781250663691529512178 ) * liVS[61] +
                   li_Real_s(       0.038483796508796741520086 ) * liVS[62] +
                   li_Real_s(      -0.005063657435894874936366 ) * liVS[63];
        liVC[45] = li_Real_s(       0.004081790124243589445996 ) * liVS[ 0] +
                   li_Real_s(      -0.026176697543412896784787 ) * liVS[ 1] +
                   li_Real_s(       0.071875000059907323546327 ) * liVS[ 2] +
                   li_Real_s(      -0.109587191490285507100566 ) * liVS[ 3] +
                   li_Real_s(       0.100270061884454642076037 ) * liVS[ 4] +
                   li_Real_s(      -0.055104166768501938733493 ) * liVS[ 5] +
                   li_Real_s(       0.016859567935937014304670 ) * liVS[ 6] +
                   li_Real_s(      -0.002218364202344669244837 ) * liVS[ 7] +
                   li_Real_s(      -0.026176697565292283975680 ) * liVS[ 8] +
                   li_Real_s(       0.167872299669256874743439 ) * liVS[ 9] +
                   li_Real_s(      -0.460937500987813497665968 ) * liVS[10] +
                   li_Real_s(       0.702787424653326997514569 ) * liVS[11] +
                   li_Real_s(      -0.643036267349164192275168 ) * liVS[12] +
                   li_Real_s(       0.353385417837228488213697 ) * liVS[13] +
                   li_Real_s(      -0.108121142359371313901306 ) * liVS[14] +
                   li_Real_s(       0.014226466101828760812964 ) * liVS[15] +
                   li_Real_s(       0.071875000158392099436355 ) * liVS[16] +
                   li_Real_s(      -0.460937501231935553569485 ) * liVS[17] +
                   li_Real_s(       1.265625004008587017878540 ) * liVS[18] +
                   li_Real_s(      -1.929687507043893468505757 ) * liVS[19] +
                   li_Real_s(       1.765625007213506014736026 ) * liVS[20] +
                   li_Real_s(      -0.970312504310568701981765 ) * liVS[21] +
                   li_Real_s(       0.296875001395175863727616 ) * liVS[22] +
                   li_Real_s(      -0.039062500189262827632319 ) * liVS[23] +
                   li_Real_s(      -0.109587191672237072026519 ) * liVS[24] +
                   li_Real_s(       0.702787425223082018455045 ) * liVS[25] +
                   li_Real_s(      -1.929687507575989391384041 ) * liVS[26] +
                   li_Real_s(       2.942177867996755402657527 ) * liVS[27] +
                   li_Real_s(      -2.692033192189889678758163 ) * liVS[28] +
                   li_Real_s(       1.479427091124125182020066 ) * liVS[29] +
                   li_Real_s(      -0.452642749416873102497050 ) * liVS[30] +
                   li_Real_s(       0.059558256511029986079997 ) * liVS[31] +
                   li_Real_s(       0.100270062065867193723534 ) * liVS[32] +
                   li_Real_s(      -0.643036267959011476058606 ) * liVS[33] +
                   li_Real_s(       1.765625007928011580560224 ) * liVS[34] +
                   li_Real_s(      -2.692033192526376073061556 ) * liVS[35] +
                   li_Real_s(       2.463155877714524155663867 ) * liVS[36] +
                   li_Real_s(      -1.353645841270910255005333 ) * liVS[37] +
                   li_Real_s(       0.414158953155060571837964 ) * liVS[38] +
                   li_Real_s(      -0.054494599107165253570884 ) * liVS[39] +
                   li_Real_s(      -0.055104166872935456122917 ) * liVS[40] +
                   li_Real_s(       0.353385418199214429435528 ) * liVS[41] +
                   li_Real_s(      -0.970312504768575778690831 ) * liVS[42] +
                   li_Real_s(       1.479427091401964267092239 ) * liVS[43] +
                   li_Real_s(      -1.353645841354618184482206 ) * liVS[44] +
                   li_Real_s(       0.743906254687345169784862 ) * liVS[45] +
                   li_Real_s(      -0.227604168159514119729181 ) * liVS[46] +
                   li_Real_s(       0.029947916867119409034537 ) * liVS[47] +
                   li_Real_s(       0.016859567969017330568704 ) * liVS[48] +
                   li_Real_s(      -0.108121142476620607508764 ) * liVS[49] +
                   li_Real_s(       0.296875001551399453347813 ) * liVS[50] +
                   li_Real_s(      -0.452642749525291820944517 ) * liVS[51] +
                   li_Real_s(       0.414158953202330426535127 ) * liVS[52] +
                   li_Real_s(      -0.227604168171723242330984 ) * liVS[53] +
                   li_Real_s(       0.069637346156869561752956 ) * liVS[54] +
                   li_Real_s(      -0.009162808705981184687062 ) * liVS[55] +
                   li_Real_s(      -0.002218364206858280951451 ) * liVS[56] +
                   li_Real_s(       0.014226466118183650477746 ) * liVS[57] +
                   li_Real_s(      -0.039062500212172696079094 ) * liVS[58] +
                   li_Real_s(       0.059558256528777220317750 ) * liVS[59] +
                   li_Real_s(      -0.054494599116627434609583 ) * liVS[60] +
                   li_Real_s(       0.029947916870565936819926 ) * liVS[61] +
                   li_Real_s(      -0.009162808706550479298514 ) * liVS[62] +
                   li_Real_s(       0.001205632724681993117599 ) * liVS[63];
        liVC[46] = li_Real_s(      -0.000354938272646376784536 ) * liVS[ 0] +
                   li_Real_s(       0.002276234575383520586378 ) * liVS[ 1] +
                   li_Real_s(      -0.006250000023296904683168 ) * liVS[ 2] +
                   li_Real_s(       0.009529321027650944353482 ) * liVS[ 3] +
                   li_Real_s(      -0.008719135842950406778407 ) * liVS[ 4] +
                   li_Real_s(       0.004791666690727156008300 ) * liVS[ 5] +
                   li_Real_s(      -0.001466049390488818882083 ) * liVS[ 6] +
                   li_Real_s(       0.000192901235621274758092 ) * liVS[ 7] +
                   li_Real_s(       0.002395833343027529982461 ) * liVS[ 8] +
                   li_Real_s(      -0.015364583402584933505275 ) * liVS[ 9] +
                   li_Real_s(       0.042187500211978257169676 ) * liVS[10] +
                   li_Real_s(      -0.064322917023614256670783 ) * liVS[11] +
                   li_Real_s(       0.058854167021844573659539 ) * liVS[12] +
                   li_Real_s(      -0.032343750208391433553778 ) * liVS[13] +
                   li_Real_s(       0.009895833400089801279442 ) * liVS[14] +
                   li_Real_s(      -0.001302083342349547034900 ) * liVS[15] +
                   li_Real_s(      -0.006921296329520959522696 ) * liVS[16] +
                   li_Real_s(       0.044386574310453785763286 ) * liVS[17] +
                   li_Real_s(      -0.121875000717080878009924 ) * liVS[18] +
                   li_Real_s(       0.185821760454857931588180 ) * liVS[19] +
                   li_Real_s(      -0.170023149327395040053545 ) * liVS[20] +
                   li_Real_s(       0.093437500687094970919233 ) * liVS[21] +
                   li_Real_s(      -0.028587963181991526973036 ) * liVS[22] +
                   li_Real_s(       0.003761574103581744044078 ) * liVS[23] +
                   li_Real_s(       0.011091821045923433430858 ) * liVS[24] +
                   li_Real_s(      -0.071132330660280754242031 ) * liVS[25] +
                   li_Real_s(       0.195312501246989927983222 ) * liVS[26] +
                   li_Real_s(      -0.297791282930733092193520 ) * liVS[27] +
                   li_Real_s(       0.272472995854292021622456 ) * liVS[28] +
                   li_Real_s(      -0.149739584509131329337350 ) * liVS[29] +
                   li_Real_s(       0.045814043583397051628836 ) * liVS[30] +
                   li_Real_s(      -0.006028163630457562469078 ) * liVS[31] +
                   li_Real_s(      -0.010648148206761853806768 ) * liVS[32] +
                   li_Real_s(       0.068287037451894905260019 ) * liVS[33] +
                   li_Real_s(      -0.187500001246365760598778 ) * liVS[34] +
                   li_Real_s(       0.285879631686053126138347 ) * liVS[35] +
                   li_Real_s(      -0.261574076083165185657720 ) * liVS[36] +
                   li_Real_s(       0.143750001161308454955545 ) * liVS[37] +
                   li_Real_s(      -0.043981481849346810986123 ) * liVS[38] +
                   li_Real_s(       0.005787037086383173267734 ) * liVS[39] +
                   li_Real_s(       0.006122685219623530405997 ) * liVS[40] +
                   li_Real_s(      -0.039265046539625300670018 ) * liVS[41] +
                   li_Real_s(       0.107812500728746318845452 ) * liVS[42] +
                   li_Real_s(      -0.164380788235331837254805 ) * liVS[43] +
                   li_Real_s(       0.150405093759580832379896 ) * liVS[44] +
                   li_Real_s(      -0.082656250672637168741552 ) * liVS[45] +
                   li_Real_s(       0.025289352064392749130928 ) * liVS[46] +
                   li_Real_s(      -0.003327546324749144045219 ) * liVS[47] +
                   li_Real_s(      -0.001952160504864478007825 ) * liVS[48] +
                   li_Real_s(       0.012519290201363410947044 ) * liVS[49] +
                   li_Real_s(      -0.034375000232793578858193 ) * liVS[50] +
                   li_Real_s(       0.052411265813903012977804 ) * liVS[51] +
                   li_Real_s(      -0.047955247284484914249703 ) * liVS[52] +
                   li_Real_s(       0.026354166879957542646284 ) * liVS[53] +
                   li_Real_s(      -0.008063271672190716610196 ) * liVS[54] +
                   li_Real_s(       0.001060956799109724624230 ) * liVS[55] +
                   li_Real_s(       0.000266203705199252738112 ) * liVS[56] +
                   li_Real_s(      -0.001707175936477305436267 ) * liVS[57] +
                   li_Real_s(       0.004687500031478691875364 ) * liVS[58] +
                   li_Real_s(      -0.007146990792270310755008 ) * liVS[59] +
                   li_Real_s(       0.006539351901814094425447 ) * liVS[60] +
                   li_Real_s(      -0.003593750028677233920860 ) * liVS[61] +
                   li_Real_s(       0.001099537046062724204853 ) * liVS[62] +
                   li_Real_s(      -0.000144675927129897519130 ) * liVS[63];
        liVC[47] = li_Real_s(       0.000012676366924407184156 ) * liVS[ 0] +
                   li_Real_s(      -0.000081294092265254015023 ) * liVS[ 1] +
                   li_Real_s(       0.000223214287351491784406 ) * liVS[ 2] +
                   li_Real_s(      -0.000340332895101973925689 ) * liVS[ 3] +
                   li_Real_s(       0.000311397709855369513032 ) * liVS[ 4] +
                   li_Real_s(      -0.000171130953902404354949 ) * liVS[ 5] +
                   li_Real_s(       0.000052358907009462289997 ) * liVS[ 6] +
                   li_Real_s(      -0.000006889329871103246419 ) * liVS[ 7] +
                   li_Real_s(      -0.000088734568547218020673 ) * liVS[ 8] +
                   li_Real_s(       0.000569058646410702416094 ) * liVS[ 9] +
                   li_Real_s(      -0.001562500013095025584720 ) * liVS[10] +
                   li_Real_s(       0.002382330268314728285994 ) * liVS[11] +
                   li_Real_s(      -0.002179783971430060331848 ) * liVS[12] +
                   li_Real_s(       0.001197916678681772931908 ) * liVS[13] +
                   li_Real_s(      -0.000366512349488999416777 ) * liVS[14] +
                   li_Real_s(       0.000048225309154099720022 ) * liVS[15] +
                   li_Real_s(       0.000266203705751909616706 ) * liVS[16] +
                   li_Real_s(      -0.001707175940028945754917 ) * liVS[17] +
                   li_Real_s(       0.004687500041587483283478 ) * liVS[18] +
                   li_Real_s(      -0.007146990808505538327056 ) * liVS[19] +
                   li_Real_s(       0.006539351917532823695034 ) * liVS[20] +
                   li_Real_s(      -0.003593750037802232854406 ) * liVS[21] +
                   li_Real_s(       0.001099537048995790754907 ) * liVS[22] +
                   li_Real_s(      -0.000144675927531290847428 ) * liVS[23] +
                   li_Real_s(      -0.000443672842949398640866 ) * liVS[24] +
                   li_Real_s(       0.002845293233606677012204 ) * liVS[25] +
                   li_Real_s(      -0.007812500069870748076184 ) * liVS[26] +
                   li_Real_s(       0.011911651348137201888999 ) * liVS[27] +
                   li_Real_s(      -0.010898919862872715955127 ) * liVS[28] +
                   li_Real_s(       0.005989583396365751565005 ) * liVS[29] +
                   li_Real_s(      -0.001832561748294010001148 ) * liVS[30] +
                   li_Real_s(       0.000241126545877252425722 ) * liVS[31] +
                   li_Real_s(       0.000443672842882688114874 ) * liVS[32] +
                   li_Real_s(      -0.002845293233146750644380 ) * liVS[33] +
                   li_Real_s(       0.007812500068402863234729 ) * liVS[34] +
                   li_Real_s(      -0.011911651345500215773421 ) * liVS[35] +
                   li_Real_s(       0.010898919860073937385114 ) * liVS[36] +
                   li_Real_s(      -0.005989583394624156849340 ) * liVS[37] +
                   li_Real_s(       0.001832561747705925515525 ) * liVS[38] +
                   li_Real_s(      -0.000241126545794291199942 ) * liVS[39] +
                   li_Real_s(      -0.000266203705655242151007 ) * liVS[40] +
                   li_Real_s(       0.001707175939371599615590 ) * liVS[41] +
                   li_Real_s(      -0.004687500039461705431121 ) * liVS[42] +
                   li_Real_s(       0.007146990804609325981245 ) * liVS[43] +
                   li_Real_s(      -0.006539351913322474611767 ) * liVS[44] +
                   li_Real_s(       0.003593750035145418849497 ) * liVS[45] +
                   li_Real_s(      -0.001099537048089650140975 ) * liVS[46] +
                   li_Real_s(       0.000144675927402726153814 ) * liVS[47] +
                   li_Real_s(       0.000088734568519977224249 ) * liVS[48] +
                   li_Real_s(      -0.000569058646237648678955 ) * liVS[49] +
                   li_Real_s(       0.001562500012493926553059 ) * liVS[50] +
                   li_Real_s(      -0.002382330267101879020530 ) * liVS[51] +
                   li_Real_s(       0.002179783970014313038144 ) * liVS[52] +
                   li_Real_s(      -0.001197916677737906852863 ) * liVS[53] +
                   li_Real_s(       0.000366512349154918780358 ) * liVS[54] +
                   li_Real_s(      -0.000048225309105700067680 ) * liVS[55] +
                   li_Real_s(      -0.000012676366926311910532 ) * liVS[56] +
                   li_Real_s(       0.000081294092284459518096 ) * liVS[57] +
                   li_Real_s(      -0.000223214287394341188986 ) * liVS[58] +
                   li_Real_s(       0.000340332895127425788528 ) * liVS[59] +
                   li_Real_s(      -0.000311397709832346263059 ) * liVS[60] +
                   li_Real_s(       0.000171130953863564438423 ) * liVS[61] +
                   li_Real_s(      -0.000052358906990368188697 ) * liVS[62] +
                   li_Real_s(       0.000006889329867915475192 ) * liVS[63];
        liVC[48] = li_Real_s(       0.005555555560625229614968 ) * liVS[ 0] +
                   li_Real_s(      -0.037500000033645321251274 ) * liVS[ 1] +
                   li_Real_s(       0.108333333429498010480607 ) * liVS[ 2] +
                   li_Real_s(      -0.173611111263739598120637 ) * liVS[ 3] +
                   li_Real_s(       0.166666666811666530634284 ) * liVS[ 4] +
                   li_Real_s(      -0.095833333415806812305426 ) * liVS[ 5] +
                   li_Real_s(       0.030555555581589465691250 ) * liVS[ 6] +
                   li_Real_s(      -0.004166666670187461375685 ) * liVS[ 7] +
                   li_Real_s(      -0.000000000031341356461395 ) * liVS[ 8] +
                   li_Real_s(       0.000000000208793849713301 ) * liVS[ 9] +
                   li_Real_s(      -0.000000000598031167566779 ) * liVS[10] +
                   li_Real_s(       0.000000000950062189287298 ) * liVS[11] +
                   li_Real_s(      -0.000000000902748569522220 ) * liVS[12] +
                   li_Real_s(       0.000000000513345447032465 ) * liVS[13] +
                   li_Real_s(      -0.000000000161975311176822 ) * liVS[14] +
                   li_Real_s(       0.000000000021894918694152 ) * liVS[15] +
                   li_Real_s(       0.000000000083230956379422 ) * liVS[16] +
                   li_Real_s(      -0.000000000556562296233379 ) * liVS[17] +
                   li_Real_s(       0.000000001597440218776323 ) * liVS[18] +
                   li_Real_s(      -0.000000002540044425928724 ) * liVS[19] +
                   li_Real_s(       0.000000002413828992644703 ) * liVS[20] +
                   li_Real_s(      -0.000000001372159113173877 ) * liVS[21] +
                   li_Real_s(       0.000000000432724864422623 ) * liVS[22] +
                   li_Real_s(      -0.000000000058459196887090 ) * liVS[23] +
                   li_Real_s(      -0.000000000122622263437733 ) * liVS[24] +
                   li_Real_s(       0.000000000823208652773864 ) * liVS[25] +
                   li_Real_s(      -0.000000002368168954690819 ) * liVS[26] +
                   li_Real_s(       0.000000003769549734394686 ) * liVS[27] +
                   li_Real_s(      -0.000000003583076089012387 ) * liVS[28] +
                   li_Real_s(       0.000000002036302222537762 ) * liVS[29] +
                   li_Real_s(      -0.000000000641855783527492 ) * liVS[30] +
                   li_Real_s(       0.000000000086662480962116 ) * liVS[31] +
                   li_Real_s(       0.000000000108142624925133 ) * liVS[32] +
                   li_Real_s(      -0.000000000729072904343114 ) * liVS[33] +
                   li_Real_s(       0.000000002102777755158136 ) * liVS[34] +
                   li_Real_s(      -0.000000003351533596127187 ) * liVS[35] +
                   li_Real_s(       0.000000003187161605991505 ) * liVS[36] +
                   li_Real_s(      -0.000000001811137999146364 ) * liVS[37] +
                   li_Real_s(       0.000000000570678335562323 ) * liVS[38] +
                   li_Real_s(      -0.000000000077015822020431 ) * liVS[39] +
                   li_Real_s(      -0.000000000057139112685303 ) * liVS[40] +
                   li_Real_s(       0.000000000386950516183424 ) * liVS[41] +
                   li_Real_s(      -0.000000001119240783862795 ) * liVS[42] +
                   li_Real_s(       0.000000001786753025567005 ) * liVS[43] +
                   li_Real_s(      -0.000000001700279419416651 ) * liVS[44] +
                   li_Real_s(       0.000000000966303783680164 ) * liVS[45] +
                   li_Real_s(      -0.000000000304416337437324 ) * liVS[46] +
                   li_Real_s(       0.000000000041068327971479 ) * liVS[47] +
                   li_Real_s(       0.000000000016783914497728 ) * liVS[48] +
                   li_Real_s(      -0.000000000114187673236994 ) * liVS[49] +
                   li_Real_s(       0.000000000331293776810085 ) * liVS[50] +
                   li_Real_s(      -0.000000000529821750697197 ) * liVS[51] +
                   li_Real_s(       0.000000000504615471709651 ) * liVS[52] +
                   li_Real_s(      -0.000000000286860726077252 ) * liVS[53] +
                   li_Real_s(       0.000000000090365277369835 ) * liVS[54] +
                   li_Real_s(      -0.000000000012188290375856 ) * liVS[55] +
                   li_Real_s(      -0.000000000002121145835574 ) * liVS[56] +
                   li_Real_s(       0.000000000014496562680919 ) * liVS[57] +
                   li_Real_s(      -0.000000000042187766232452 ) * liVS[58] +
                   li_Real_s(       0.000000000067593192320045 ) * liVS[59] +
                   li_Real_s(      -0.000000000064438413707393 ) * liVS[60] +
                   li_Real_s(       0.000000000036644767088283 ) * liVS[61] +
                   li_Real_s(      -0.000000000011544015833222 ) * liVS[62] +
                   li_Real_s(       0.000000000001556819519396 ) * liVS[63];
        liVC[49] = li_Real_s(      -0.014404761898050733037735 ) * liVS[ 0] +
                   li_Real_s(       0.097232142819622002782864 ) * liVS[ 1] +
                   li_Real_s(      -0.280892857059469591707312 ) * liVS[ 2] +
                   li_Real_s(       0.450148809428364327089866 ) * liVS[ 3] +
                   li_Real_s(      -0.432142857082866815865430 ) * liVS[ 4] +
                   li_Real_s(       0.248482142837158748172044 ) * liVS[ 5] +
                   li_Real_s(      -0.079226190473214747056474 ) * liVS[ 6] +
                   li_Real_s(       0.010803571428458169645381 ) * liVS[ 7] +
                   li_Real_s(       0.038888888865275972328561 ) * liVS[ 8] +
                   li_Real_s(      -0.262499999894393598598441 ) * liVS[ 9] +
                   li_Real_s(       0.758333333189271763252748 ) * liVS[10] +
                   li_Real_s(      -1.215277777782534540662596 ) * liVS[11] +
                   li_Real_s(       1.166666666859635048680843 ) * liVS[12] +
                   li_Real_s(      -0.670833333527427466691506 ) * liVS[13] +
                   li_Real_s(       0.213888888969268775586130 ) * liVS[14] +
                   li_Real_s(      -0.029166666679095953895740 ) * liVS[15] +
                   li_Real_s(      -0.058333333312040824836231 ) * liVS[16] +
                   li_Real_s(       0.393750000011787004705610 ) * liVS[17] +
                   li_Real_s(      -1.137500000474837014152740 ) * liVS[18] +
                   li_Real_s(       1.822916668049823440966861 ) * liVS[19] +
                   li_Real_s(      -1.750000001809369543082084 ) * liVS[20] +
                   li_Real_s(       1.006250001244665570609982 ) * liVS[21] +
                   li_Real_s(      -0.320833333773019757639844 ) * liVS[22] +
                   li_Real_s(       0.043750000062990901383841 ) * liVS[23] +
                   li_Real_s(       0.064814814832274247891064 ) * liVS[24] +
                   li_Real_s(      -0.437500000365467323071300 ) * liVS[25] +
                   li_Real_s(       1.263888890615923132898502 ) * liVS[26] +
                   li_Real_s(      -2.025462966649275831088062 ) * liVS[27] +
                   li_Real_s(       1.944444448673531722349139 ) * liVS[28] +
                   li_Real_s(      -1.118055558268958860068665 ) * liVS[29] +
                   li_Real_s(       0.356481482401948768234945 ) * liVS[30] +
                   li_Real_s(      -0.048611111239975475506458 ) * liVS[31] +
                   li_Real_s(      -0.048611111161257047896811 ) * liVS[32] +
                   li_Real_s(       0.328125000573676328663453 ) * liVS[33] +
                   li_Real_s(      -0.947916668938157291890434 ) * liVS[34] +
                   li_Real_s(       1.519097226680505174911673 ) * liVS[35] +
                   li_Real_s(      -1.458333338211462493205772 ) * liVS[36] +
                   li_Real_s(       0.838541669706514780635587 ) * liVS[37] +
                   li_Real_s(      -0.267361112123226218884042 ) * liVS[38] +
                   li_Real_s(       0.036458333473406989710952 ) * liVS[39] +
                   li_Real_s(       0.023333333375125775432934 ) * liVS[40] +
                   li_Real_s(      -0.157500000415599528125199 ) * liVS[41] +
                   li_Real_s(       0.455000001538358955777142 ) * liVS[42] +
                   li_Real_s(      -0.729166669571096770496865 ) * liVS[43] +
                   li_Real_s(       0.700000003101284073991906 ) * liVS[44] +
                   li_Real_s(      -0.402500001901542669990874 ) * liVS[45] +
                   li_Real_s(       0.128333333959503470556029 ) * liVS[46] +
                   li_Real_s(      -0.017500000086034108587318 ) * liVS[47] +
                   li_Real_s(      -0.006481481497545349412803 ) * liVS[48] +
                   li_Real_s(       0.043750000150423462574878 ) * liVS[49] +
                   li_Real_s(      -0.126388889426850559161153 ) * liVS[50] +
                   li_Real_s(       0.202546297290010501868096 ) * liVS[51] +
                   li_Real_s(      -0.194444445490042916446782 ) * liVS[52] +
                   li_Real_s(       0.111805556190114452341788 ) * liVS[53] +
                   li_Real_s(      -0.035648148355574948986657 ) * liVS[54] +
                   li_Real_s(       0.004861111139465523756087 ) * liVS[55] +
                   li_Real_s(       0.000793650796067968400394 ) * liVS[56] +
                   li_Real_s(      -0.005357142879110855793190 ) * liVS[57] +
                   li_Real_s(       0.015476190553251778503352 ) * liVS[58] +
                   li_Real_s(      -0.024801587442071393319054 ) * liVS[59] +
                   li_Real_s(       0.023809523955971079178795 ) * liVS[60] +
                   li_Real_s(      -0.013690476278746649060891 ) * liVS[61] +
                   li_Real_s(       0.004365079393784693229108 ) * liVS[62] +
                   li_Real_s(      -0.000595238099147077370787 ) * liVS[63];
        liVC[50] = li_Real_s(       0.014475308612236226224468 ) * liVS[ 0] +
                   li_Real_s(      -0.097708333153270876536567 ) * liVS[ 1] +
                   li_Real_s(       0.282268518060947104508784 ) * liVS[ 2] +
                   li_Real_s(      -0.452353394422889465431581 ) * liVS[ 3] +
                   li_Real_s(       0.434259258724592367428841 ) * liVS[ 4] +
                   li_Real_s(      -0.249699073803025584084025 ) * liVS[ 5] +
                   li_Real_s(       0.079614197453016982697349 ) * liVS[ 6] +
                   li_Real_s(      -0.010856481471610557321128 ) * liVS[ 7] +
                   li_Real_s(      -0.061944444301423651211280 ) * liVS[ 8] +
                   li_Real_s(       0.418124999173804912100394 ) * liVS[ 9] +
                   li_Real_s(      -1.207916664679933749226848 ) * liVS[10] +
                   li_Real_s(       1.935763886287902746374812 ) * liVS[11] +
                   li_Real_s(      -1.858333331307619573635748 ) * liVS[12] +
                   li_Real_s(       1.068541665713139243720775 ) * liVS[13] +
                   li_Real_s(      -0.340694444188130551864901 ) * liVS[14] +
                   li_Real_s(       0.046458333302260790276250 ) * liVS[15] +
                   li_Real_s(       0.122083333045443165332244 ) * liVS[16] +
                   li_Real_s(      -0.824062498451672453825267 ) * liVS[17] +
                   li_Real_s(       2.380624996615860133886144 ) * liVS[18] +
                   li_Real_s(      -3.815104162785000418978143 ) * liVS[19] +
                   li_Real_s(       3.662499997488128933298412 ) * liVS[20] +
                   li_Real_s(      -2.105937499083414099487754 ) * liVS[21] +
                   li_Real_s(       0.671458333154790865116013 ) * liVS[22] +
                   li_Real_s(      -0.091562499984136458408557 ) * liVS[23] +
                   li_Real_s(      -0.146450616965413171755017 ) * liVS[24] +
                   li_Real_s(       0.988541665134555547211903 ) * liVS[25] +
                   li_Real_s(      -2.855787034257553980154398 ) * liVS[26] +
                   li_Real_s(       4.576581787944891566155547 ) * liVS[27] +
                   li_Real_s(      -4.393518518179516441080068 ) * liVS[28] +
                   li_Real_s(       2.526273148692343539778449 ) * liVS[29] +
                   li_Real_s(      -0.805478395386727274996019 ) * liVS[30] +
                   li_Real_s(       0.109837963017422504674592 ) * liVS[31] +
                   li_Real_s(       0.113888888681118682910665 ) * liVS[32] +
                   li_Real_s(      -0.768749999174872744589493 ) * liVS[33] +
                   li_Real_s(       2.220833332443664431821162 ) * liVS[34] +
                   li_Real_s(      -3.559027778373651251797583 ) * liVS[35] +
                   li_Real_s(       3.416666668763476355508146 ) * liVS[36] +
                   li_Real_s(      -1.964583335163705424264435 ) * liVS[37] +
                   li_Real_s(       0.626388889593805009425864 ) * liVS[38] +
                   li_Real_s(      -0.085416666769835503103536 ) * liVS[39] +
                   li_Real_s(      -0.055833333254826644775903 ) * liVS[40] +
                   li_Real_s(       0.376874999792727871650300 ) * liVS[41] +
                   li_Real_s(      -1.088750000215061186636945 ) * liVS[42] +
                   li_Real_s(       1.744791668091729253120548 ) * liVS[43] +
                   li_Real_s(      -1.675000002167238655204073 ) * liVS[44] +
                   li_Real_s(       0.963125001564241522977738 ) * liVS[45] +
                   li_Real_s(      -0.307083333890290410828072 ) * liVS[46] +
                   li_Real_s(       0.041875000078718374596498 ) * liVS[47] +
                   li_Real_s(       0.015725308626392875410716 ) * liVS[48] +
                   li_Real_s(      -0.106145833328253280658515 ) * liVS[49] +
                   li_Real_s(       0.306643518750683519158429 ) * liVS[50] +
                   li_Real_s(      -0.491415895755851206416764 ) * liVS[51] +
                   li_Real_s(       0.471759260159539728363143 ) * liVS[52] +
                   li_Real_s(      -0.271261574682206574316012 ) * liVS[53] +
                   li_Real_s(       0.086489197740280854276307 ) * liVS[54] +
                   li_Real_s(      -0.011793981510585638261546 ) * liVS[55] +
                   li_Real_s(      -0.001944444443206627681775 ) * liVS[56] +
                   li_Real_s(       0.013125000004950239385071 ) * liVS[57] +
                   li_Real_s(      -0.037916666713147417766550 ) * liVS[58] +
                   li_Real_s(       0.060763889004726401310563 ) * liVS[59] +
                   li_Real_s(      -0.058333333474076543012643 ) * liVS[60] +
                   li_Real_s(       0.033541666758708538198519 ) * liVS[61] +
                   li_Real_s(      -0.010694444475570968888789 ) * liVS[62] +
                   li_Real_s(       0.001458333337615913549712 ) * liVS[63];
        liVC[51] = li_Real_s(      -0.007461419732555540917929 ) * liVS[ 0] +
                   li_Real_s(       0.050364583211120139694117 ) * liVS[ 1] +
                   li_Real_s(      -0.145497684881120914468511 ) * liVS[ 2] +
                   li_Real_s(       0.233169366869966143696047 ) * liVS[ 3] +
                   li_Real_s(      -0.223842592255599326378501 ) * liVS[ 4] +
                   li_Real_s(       0.128709490574477891655647 ) * liVS[ 5] +
                   li_Real_s(      -0.041037808595253799648361 ) * liVS[ 6] +
                   li_Real_s(       0.005596064808958689518192 ) * liVS[ 7] +
                   li_Real_s(       0.039382715956428282488844 ) * liVS[ 8] +
                   li_Real_s(      -0.265833332817610401566100 ) * liVS[ 9] +
                   li_Real_s(       0.767962961788526410344957 ) * liVS[10] +
                   li_Real_s(      -1.230709875115879370355287 ) * liVS[11] +
                   li_Real_s(       1.181481480474605438857338 ) * liVS[12] +
                   li_Real_s(      -0.679351851432567599431422 ) * liVS[13] +
                   li_Real_s(       0.216604938172624772452224 ) * liVS[14] +
                   li_Real_s(      -0.029537037026127421768251 ) * liVS[15] +
                   li_Real_s(      -0.090949073903199950663634 ) * liVS[16] +
                   li_Real_s(       0.613906249165464323880315 ) * liVS[17] +
                   li_Real_s(      -1.773506942896516447305544 ) * liVS[18] +
                   li_Real_s(       2.842158563541798077523026 ) * liVS[19] +
                   li_Real_s(      -2.728472221938831321352836 ) * liVS[20] +
                   li_Real_s(       1.568871528013223048958480 ) * liVS[21] +
                   li_Real_s(      -0.500219907564276256017877 ) * liVS[22] +
                   li_Real_s(       0.068211805582338858044977 ) * liVS[23] +
                   li_Real_s(       0.120061728230859365851302 ) * liVS[24] +
                   li_Real_s(      -0.810416666058686008256018 ) * liVS[25] +
                   li_Real_s(       2.341203703255875012700926 ) * liVS[26] +
                   li_Real_s(      -3.751929013414227576106441 ) * liVS[27] +
                   li_Real_s(       3.601851854232957172996521 ) * liVS[28] +
                   li_Real_s(      -2.071064816740239233894272 ) * liVS[29] +
                   li_Real_s(       0.660339506894093597466622 ) * liVS[30] +
                   li_Real_s(      -0.090046296400627459655119 ) * liVS[31] +
                   li_Real_s(      -0.098186728311540960589809 ) * liVS[32] +
                   li_Real_s(       0.662760416570783994494320 ) * liVS[33] +
                   li_Real_s(      -1.914641204614865310773553 ) * liVS[34] +
                   li_Real_s(       3.068335265392852395649470 ) * liVS[35] +
                   li_Real_s(      -2.945601855954943282256409 ) * liVS[36] +
                   li_Real_s(       1.693721067641771327316746 ) * liVS[37] +
                   li_Real_s(      -0.540027007158383476337349 ) * liVS[38] +
                   li_Real_s(       0.073640046434325645563490 ) * liVS[39] +
                   li_Real_s(       0.049629629612089276591291 ) * liVS[40] +
                   li_Real_s(      -0.335000000143070686409175 ) * liVS[41] +
                   li_Real_s(       0.967777778838110558368157 ) * liVS[42] +
                   li_Real_s(      -1.550925928472464132568120 ) * liVS[43] +
                   li_Real_s(       1.488888891949352810684104 ) * liVS[44] +
                   li_Real_s(      -0.856111113107281829215367 ) * liVS[45] +
                   li_Real_s(       0.272962963638716638037351 ) * liVS[46] +
                   li_Real_s(      -0.037222222315453135088603 ) * liVS[47] +
                   li_Real_s(      -0.014266975309933482840563 ) * liVS[48] +
                   li_Real_s(       0.096302083423171730425238 ) * liVS[49] +
                   li_Real_s(      -0.278206018974194169146585 ) * liVS[50] +
                   li_Real_s(       0.445842979377559789355701 ) * liVS[51] +
                   li_Real_s(      -0.428009260382178347015270 ) * liVS[52] +
                   li_Real_s(       0.246105324786440649020847 ) * liVS[53] +
                   li_Real_s(      -0.078468364434779447336155 ) * liVS[54] +
                   li_Real_s(       0.010700231513913527336967 ) * liVS[55] +
                   li_Real_s(       0.001790123457592773803526 ) * liVS[56] +
                   li_Real_s(      -0.012083333349532293654605 ) * liVS[57] +
                   li_Real_s(       0.034907407479762730950767 ) * liVS[58] +
                   li_Real_s(      -0.055941358172990840458283 ) * liVS[59] +
                   li_Real_s(       0.053703703868697605372518 ) * liVS[60] +
                   li_Real_s(      -0.030879629732623030402561 ) * liVS[61] +
                   li_Real_s(       0.009845679046296629266521 ) * liVS[62] +
                   li_Real_s(      -0.001342592597202596493844 ) * liVS[63];
        liVC[52] = li_Real_s(       0.002160493822764841809203 ) * liVS[ 0] +
                   li_Real_s(      -0.014583333309301782509948 ) * liVS[ 1] +
                   li_Real_s(       0.042129629577156824105089 ) * liVS[ 2] +
                   li_Real_s(      -0.067515432039971567945713 ) * liVS[ 3] +
                   li_Real_s(       0.064814814778965612607209 ) * liVS[ 4] +
                   li_Real_s(      -0.037268518507024395169935 ) * liVS[ 5] +
                   li_Real_s(       0.011882716047691650740337 ) * liVS[ 6] +
                   li_Real_s(      -0.001620370370280510563532 ) * liVS[ 7] +
                   li_Real_s(      -0.012847222208953068545156 ) * liVS[ 8] +
                   li_Real_s(       0.086718749949846851698965 ) * liVS[ 9] +
                   li_Real_s(      -0.250520833299170242458587 ) * liVS[10] +
                   li_Real_s(       0.401475694547848038773452 ) * liVS[11] +
                   li_Real_s(      -0.385416666889929149242278 ) * liVS[12] +
                   li_Real_s(       0.221614583514228524752099 ) * liVS[13] +
                   li_Real_s(      -0.070659722290486060014558 ) * liVS[14] +
                   li_Real_s(       0.009635416676615049524912 ) * liVS[15] +
                   li_Real_s(       0.032870370366385737170845 ) * liVS[16] +
                   li_Real_s(      -0.221875000087362550527104 ) * liVS[17] +
                   li_Real_s(       0.640972222780224032412377 ) * liVS[18] +
                   li_Real_s(      -1.027199075389977123151652 ) * liVS[19] +
                   li_Real_s(       0.986111112694279978718725 ) * liVS[20] +
                   li_Real_s(      -0.567013889928779657445546 ) * liVS[21] +
                   li_Real_s(       0.180787037392591609652470 ) * liVS[22] +
                   li_Real_s(      -0.024652777827362082341267 ) * liVS[23] +
                   li_Real_s(      -0.047029321018673631016327 ) * liVS[24] +
                   li_Real_s(       0.317447917056945472236151 ) * liVS[25] +
                   li_Real_s(      -0.917071760852253614615393 ) * liVS[26] +
                   li_Real_s(       1.469666284041232895418716 ) * liVS[27] +
                   li_Real_s(      -1.410879633139494826465921 ) * liVS[28] +
                   li_Real_s(       0.811255789233390700587734 ) * liVS[29] +
                   li_Real_s(      -0.258661266162247494015958 ) * liVS[30] +
                   li_Real_s(       0.035271990841100421543164 ) * liVS[31] +
                   li_Real_s(       0.040740740794886320941259 ) * liVS[32] +
                   li_Real_s(      -0.275000000535428446024611 ) * liVS[33] +
                   li_Real_s(       0.794444446414484972684988 ) * liVS[34] +
                   li_Real_s(      -1.273148151858678689407611 ) * liVS[35] +
                   li_Real_s(       1.222222226181172377579287 ) * liVS[36] +
                   li_Real_s(      -0.702777780201011226246521 ) * liVS[37] +
                   li_Real_s(       0.224074074868443107177995 ) * liVS[38] +
                   li_Real_s(      -0.030555555663868583238241 ) * liVS[39] +
                   li_Real_s(      -0.021412037077750190050551 ) * liVS[40] +
                   li_Real_s(       0.144531250370051378428826 ) * liVS[41] +
                   li_Real_s(      -0.417534723521672146429040 ) * liVS[42] +
                   li_Real_s(       0.669126159784272367403446 ) * liVS[43] +
                   li_Real_s(      -0.642361113598520105938405 ) * liVS[44] +
                   li_Real_s(       0.369357640391594299611455 ) * liVS[45] +
                   li_Real_s(      -0.117766204191955076163723 ) * liVS[46] +
                   li_Real_s(       0.016059027843979226807258 ) * liVS[47] +
                   li_Real_s(       0.006327160508841134234359 ) * liVS[48] +
                   li_Real_s(      -0.042708333464288379677232 ) * liVS[49] +
                   li_Real_s(       0.123379630077751656358487 ) * liVS[50] +
                   li_Real_s(      -0.197723766237563158831847 ) * liVS[51] +
                   li_Real_s(       0.189814815647507240115033 ) * liVS[52] +
                   li_Real_s(      -0.109143519017210532950912 ) * liVS[53] +
                   li_Real_s(       0.034799382877070850206280 ) * liVS[54] +
                   li_Real_s(      -0.004745370392108816393062 ) * liVS[55] +
                   li_Real_s(      -0.000810185187398615447307 ) * liVS[56] +
                   li_Real_s(       0.005468750018892083730737 ) * liVS[57] +
                   li_Real_s(      -0.015798611174778487420411 ) * liVS[58] +
                   li_Real_s(       0.025318287150224882964267 ) * liVS[59] +
                   li_Real_s(      -0.024305555671631118297427 ) * liVS[60] +
                   li_Real_s(       0.013975694513543034797065 ) * liVS[61] +
                   li_Real_s(      -0.004456018540726608412328 ) * liVS[62] +
                   li_Real_s(       0.000607638891875927900088 ) * liVS[63];
        liVC[53] = li_Real_s(      -0.000354938271682592176859 ) * liVS[ 0] +
                   li_Real_s(       0.002395833334752767163067 ) * liVS[ 1] +
                   li_Real_s(      -0.006921296303176692477077 ) * liVS[ 2] +
                   li_Real_s(       0.011091821002854163058515 ) * liVS[ 3] +
                   li_Real_s(      -0.010648148166083268306714 ) * liVS[ 4] +
                   li_Real_s(       0.006122685196912724947538 ) * liVS[ 5] +
                   li_Real_s(      -0.001952160497842171660299 ) * liVS[ 6] +
                   li_Real_s(       0.000266203704264923735057 ) * liVS[ 7] +
                   li_Real_s(       0.002276234571402802053797 ) * liVS[ 8] +
                   li_Real_s(      -0.015364583363523169265430 ) * liVS[ 9] +
                   li_Real_s(       0.044386574179864787814687 ) * liVS[10] +
                   li_Real_s(      -0.071132330442795155267532 ) * liVS[11] +
                   li_Real_s(       0.068287037245171197663574 ) * liVS[12] +
                   li_Real_s(      -0.039265046423952475185271 ) * liVS[13] +
                   li_Real_s(       0.012519290165531285485834 ) * liVS[14] +
                   li_Real_s(      -0.001707175931699278503828 ) * liVS[15] +
                   li_Real_s(      -0.006250000016176926465050 ) * liVS[16] +
                   li_Real_s(       0.042187500128688493195028 ) * liVS[17] +
                   li_Real_s(      -0.121875000423684876071917 ) * liVS[18] +
                   li_Real_s(       0.195312500749458217708110 ) * liVS[19] +
                   li_Real_s(      -0.187500000770926378290682 ) * liVS[20] +
                   li_Real_s(       0.107812500462388161093230 ) * liVS[21] +
                   li_Real_s(      -0.034375000150194123804681 ) * liVS[22] +
                   li_Real_s(       0.004687500020447432635962 ) * liVS[23] +
                   li_Real_s(       0.009529321019878800314018 ) * liVS[24] +
                   li_Real_s(      -0.064322916915478950405927 ) * liVS[25] +
                   li_Real_s(       0.185821760057548718281240 ) * liVS[26] +
                   li_Real_s(      -0.297791282247609978561798 ) * liVS[27] +
                   li_Real_s(       0.285879631030745040565222 ) * liVS[28] +
                   li_Real_s(      -0.164380787867923733314157 ) * liVS[29] +
                   li_Real_s(       0.052411265699860321021220 ) * liVS[30] +
                   li_Real_s(      -0.007146990777019899578060 ) * liVS[31] +
                   li_Real_s(      -0.008719135837252589560364 ) * liVS[32] +
                   li_Real_s(       0.058854166931067361523411 ) * liVS[33] +
                   li_Real_s(      -0.170023148984264460548133 ) * liVS[34] +
                   li_Real_s(       0.272472995258786709893428 ) * liVS[35] +
                   li_Real_s(      -0.261574075510336068184358 ) * liVS[36] +
                   li_Real_s(       0.150405093438165909924820 ) * liVS[37] +
                   li_Real_s(      -0.047955247184606704380361 ) * liVS[38] +
                   li_Real_s(       0.006539351888439820514876 ) * liVS[39] +
                   li_Real_s(       0.004791666688049235622859 ) * liVS[40] +
                   li_Real_s(      -0.032343750160991169284586 ) * liVS[41] +
                   li_Real_s(       0.093437500504347528740290 ) * liVS[42] +
                   li_Real_s(      -0.149739584189920477141911 ) * liVS[43] +
                   li_Real_s(       0.143750000853700848901795 ) * liVS[44] +
                   li_Real_s(      -0.082656250499960420841106 ) * liVS[45] +
                   li_Real_s(       0.026354166826252141409004 ) * liVS[46] +
                   li_Real_s(      -0.003593750021477735978603 ) * liVS[47] +
                   li_Real_s(      -0.001466049389786470980024 ) * liVS[48] +
                   li_Real_s(       0.009895833386252942537320 ) * liVS[49] +
                   li_Real_s(      -0.028587963127703869270491 ) * liVS[50] +
                   li_Real_s(       0.045814043488085459676640 ) * liVS[51] +
                   li_Real_s(      -0.043981481757424334411155 ) * liVS[52] +
                   li_Real_s(       0.025289352012812051084811 ) * liVS[53] +
                   li_Real_s(      -0.008063271656146758781247 ) * liVS[54] +
                   li_Real_s(       0.001099537043910983613593 ) * liVS[55] +
                   li_Real_s(       0.000192901235546938387699 ) * liVS[56] +
                   li_Real_s(      -0.001302083340635740420915 ) * liVS[57] +
                   li_Real_s(       0.003761574096710454048753 ) * liVS[58] +
                   li_Real_s(      -0.006028163618321202843475 ) * liVS[59] +
                   li_Real_s(       0.005787037074668398284327 ) * liVS[60] +
                   li_Real_s(      -0.003327546318180101860007 ) * liVS[61] +
                   li_Real_s(       0.001060956797067082527097 ) * liVS[62] +
                   li_Real_s(      -0.000144675926856019376743 ) * liVS[63];
        liVC[54] = li_Real_s(       0.000030864197648101249216 ) * liVS[ 0] +
                   li_Real_s(      -0.000208333334183491425406 ) * liVS[ 1] +
                   li_Real_s(       0.000601851854507803527006 ) * liVS[ 2] +
                   li_Real_s(      -0.000964506177402780356944 ) * liVS[ 3] +
                   li_Real_s(       0.000925925930546538447641 ) * liVS[ 4] +
                   li_Real_s(      -0.000532407410157097887859 ) * liVS[ 5] +
                   li_Real_s(       0.000169753087310217567335 ) * liVS[ 6] +
                   li_Real_s(      -0.000023148148269275942157 ) * liVS[ 7] +
                   li_Real_s(      -0.000208333334392256722123 ) * liVS[ 8] +
                   li_Real_s(       0.001406250007621641980732 ) * liVS[ 9] +
                   li_Real_s(      -0.004062500023408346580545 ) * liVS[10] +
                   li_Real_s(       0.006510416706143527104667 ) * liVS[11] +
                   li_Real_s(      -0.006250000039324147584063 ) * liVS[12] +
                   li_Real_s(       0.003593750023105205405188 ) * liVS[13] +
                   li_Real_s(      -0.001145833340749976593520 ) * liVS[14] +
                   li_Real_s(       0.000156250001004352989664 ) * liVS[15] +
                   li_Real_s(       0.000601851855435261290372 ) * liVS[16] +
                   li_Real_s(      -0.004062500025665510654249 ) * liVS[17] +
                   li_Real_s(       0.011736111189217109454508 ) * liVS[18] +
                   li_Real_s(      -0.018807870500800166063682 ) * liVS[19] +
                   li_Real_s(       0.018055555684343846389872 ) * liVS[20] +
                   li_Real_s(      -0.010381944519580016753069 ) * liVS[21] +
                   li_Real_s(       0.003310185209177322375462 ) * liVS[22] +
                   li_Real_s(      -0.000451388892127849508662 ) * liVS[23] +
                   li_Real_s(      -0.000964506179090537929532 ) * liVS[24] +
                   li_Real_s(       0.006510416711290252164690 ) * liVS[25] +
                   li_Real_s(      -0.018807870505386039411411 ) * liVS[26] +
                   li_Real_s(       0.030140818125310893382807 ) * liVS[27] +
                   li_Real_s(      -0.028935185405209694342599 ) * liVS[28] +
                   li_Real_s(       0.016637731609246857178430 ) * liVS[29] +
                   li_Real_s(      -0.005304783991265231771206 ) * liVS[30] +
                   li_Real_s(       0.000723379635103484140528 ) * liVS[31] +
                   li_Real_s(       0.000925925932204463980613 ) * liVS[32] +
                   li_Real_s(      -0.006250000044703950957681 ) * liVS[33] +
                   li_Real_s(       0.018055555690235154664425 ) * liVS[34] +
                   li_Real_s(      -0.028935185407691202397196 ) * liVS[35] +
                   li_Real_s(       0.027777777995354340639933 ) * liVS[36] +
                   li_Real_s(      -0.015972222348113570417061 ) * liVS[37] +
                   li_Real_s(       0.005092592632524419946072 ) * liVS[38] +
                   li_Real_s(      -0.000694444449809658928552 ) * liVS[39] +
                   li_Real_s(      -0.000532407411099132532595 ) * liVS[40] +
                   li_Real_s(       0.003593750026232980787633 ) * liVS[41] +
                   li_Real_s(      -0.010381944523213107109405 ) * liVS[42] +
                   li_Real_s(       0.016637731611153859512253 ) * liVS[43] +
                   li_Real_s(      -0.015972222348602412023144 ) * liVS[44] +
                   li_Real_s(       0.009184027850686985253170 ) * liVS[45] +
                   li_Real_s(      -0.002928240763807441263133 ) * liVS[46] +
                   li_Real_s(       0.000399305558648263797354 ) * liVS[47] +
                   li_Real_s(       0.000169753087605566815865 ) * liVS[48] +
                   li_Real_s(      -0.001145833341748604640095 ) * liVS[49] +
                   li_Real_s(       0.003310185210390666904479 ) * liVS[50] +
                   li_Real_s(      -0.005304783991998057030015 ) * liVS[51] +
                   li_Real_s(       0.005092592632817243003540 ) * liVS[52] +
                   li_Real_s(      -0.002928240763890651611467 ) * liVS[53] +
                   li_Real_s(       0.000933641982616644980353 ) * liVS[54] +
                   li_Real_s(      -0.000127314815792806471095 ) * liVS[55] +
                   li_Real_s(      -0.000023148148309334176664 ) * liVS[56] +
                   li_Real_s(       0.000156250001143053455707 ) * liVS[57] +
                   li_Real_s(      -0.000451388892306370768936 ) * liVS[58] +
                   li_Real_s(       0.000723379635228530923352 ) * liVS[59] +
                   li_Real_s(      -0.000694444449875736280475 ) * liVS[60] +
                   li_Real_s(       0.000399305558675219340187 ) * liVS[61] +
                   li_Real_s(      -0.000127314815797798571578 ) * liVS[62] +
                   li_Real_s(       0.000017361111242435861568 ) * liVS[63];
        liVC[55] = li_Real_s(      -0.000001102292778156974129 ) * liVS[ 0] +
                   li_Real_s(       0.000007440476253469275204 ) * liVS[ 1] +
                   li_Real_s(      -0.000021494709181047946681 ) * liVS[ 2] +
                   li_Real_s(       0.000034446649335792317755 ) * liVS[ 3] +
                   li_Real_s(      -0.000033068783367673818506 ) * liVS[ 4] +
                   li_Real_s(       0.000019014550437966614219 ) * liVS[ 5] +
                   li_Real_s(      -0.000006062610284519931501 ) * liVS[ 6] +
                   li_Real_s(       0.000000826719584171134491 ) * liVS[ 7] +
                   li_Real_s(       0.000007716049454307001522 ) * liVS[ 8] +
                   li_Real_s(      -0.000052083333826747260490 ) * liVS[ 9] +
                   li_Real_s(       0.000150462964421747022545 ) * liVS[10] +
                   li_Real_s(      -0.000241126545594823370536 ) * liVS[11] +
                   li_Real_s(       0.000231481483801414147754 ) * liVS[12] +
                   li_Real_s(      -0.000133101853191981743876 ) * liVS[13] +
                   li_Real_s(       0.000042438272030413117745 ) * liVS[14] +
                   li_Real_s(      -0.000005787037094328846901 ) * liVS[15] +
                   li_Real_s(      -0.000023148148372681182257 ) * liVS[16] +
                   li_Real_s(       0.000156250001551339978787 ) * liVS[17] +
                   li_Real_s(      -0.000451388893469464207771 ) * liVS[18] +
                   li_Real_s(       0.000723379637096151564823 ) * liVS[19] +
                   li_Real_s(      -0.000694444451682964150434 ) * liVS[20] +
                   li_Real_s(       0.000399305559723841358398 ) * liVS[21] +
                   li_Real_s(      -0.000127314816134781996684 ) * liVS[22] +
                   li_Real_s(       0.000017361111288558472507 ) * liVS[23] +
                   li_Real_s(       0.000038580247288905362979 ) * liVS[24] +
                   li_Real_s(      -0.000260416669262033943517 ) * liVS[25] +
                   li_Real_s(       0.000752314822466388852251 ) * liVS[26] +
                   li_Real_s(      -0.001205632728490920683312 ) * liVS[27] +
                   li_Real_s(       0.001157407419436920579386 ) * liVS[28] +
                   li_Real_s(      -0.000665509266169051968975 ) * liVS[29] +
                   li_Real_s(       0.000212191360208007222428 ) * liVS[30] +
                   li_Real_s(      -0.000028935185478216678237 ) * liVS[31] +
                   li_Real_s(      -0.000038580247280720070258 ) * liVS[32] +
                   li_Real_s(       0.000260416669205379338464 ) * liVS[33] +
                   li_Real_s(      -0.000752314822286924444507 ) * liVS[34] +
                   li_Real_s(       0.001205632728171590834290 ) * liVS[35] +
                   li_Real_s(      -0.001157407419100621882163 ) * liVS[36] +
                   li_Real_s(       0.000665509265960871274070 ) * liVS[37] +
                   li_Real_s(      -0.000212191360137909350693 ) * liVS[38] +
                   li_Real_s(       0.000028935185468334165271 ) * liVS[39] +
                   li_Real_s(       0.000023148148360246467541 ) * liVS[40] +
                   li_Real_s(      -0.000156250001466181210729 ) * liVS[41] +
                   li_Real_s(       0.000451388893197413433936 ) * liVS[42] +
                   li_Real_s(      -0.000723379636605280765954 ) * liVS[43] +
                   li_Real_s(       0.000694444451159199741853 ) * liVS[44] +
                   li_Real_s(      -0.000399305559396175946234 ) * liVS[45] +
                   li_Real_s(       0.000127314816023584975041 ) * liVS[46] +
                   li_Real_s(      -0.000017361111272806532824 ) * liVS[47] +
                   li_Real_s(      -0.000007716049450025866613 ) * liVS[48] +
                   li_Real_s(       0.000052083333798535310299 ) * liVS[49] +
                   li_Real_s(      -0.000150462964328237675145 ) * liVS[50] +
                   li_Real_s(       0.000241126545416591375401 ) * liVS[51] +
                   li_Real_s(      -0.000231481483601941005757 ) * liVS[52] +
                   li_Real_s(       0.000133101853062587201920 ) * liVS[53] +
                   li_Real_s(      -0.000042438271985356982273 ) * liVS[54] +
                   li_Real_s(       0.000005787037087847594734 ) * liVS[55] +
                   li_Real_s(       0.000001102292778037495050 ) * liVS[56] +
                   li_Real_s(      -0.000007440476253207789357 ) * liVS[57] +
                   li_Real_s(       0.000021494709178625784817 ) * liVS[58] +
                   li_Real_s(      -0.000034446649326846565629 ) * liVS[59] +
                   li_Real_s(       0.000033068783353631963804 ) * liVS[60] +
                   li_Real_s(      -0.000019014550426954406634 ) * liVS[61] +
                   li_Real_s(       0.000006062610280230075542 ) * liVS[62] +
                   li_Real_s(      -0.000000826719583516032433 ) * liVS[63];
        liVC[56] = li_Real_s(      -0.000198412698629489380299 ) * liVS[ 0] +
                   li_Real_s(       0.001388888890331262506231 ) * liVS[ 1] +
                   li_Real_s(      -0.004166666670796116862807 ) * liVS[ 2] +
                   li_Real_s(       0.006944444451005343452166 ) * liVS[ 3] +
                   li_Real_s(      -0.006944444450681311852003 ) * liVS[ 4] +
                   li_Real_s(       0.004166666670215458077864 ) * liVS[ 5] +
                   li_Real_s(      -0.001388888890009412103999 ) * liVS[ 6] +
                   li_Real_s(       0.000198412698564267111524 ) * liVS[ 7] +
                   li_Real_s(       0.000000000001336826352772 ) * liVS[ 8] +
                   li_Real_s(      -0.000000000008929060839562 ) * liVS[ 9] +
                   li_Real_s(       0.000000000025619763987780 ) * liVS[10] +
                   li_Real_s(      -0.000000000040746473755518 ) * liVS[11] +
                   li_Real_s(       0.000000000038744327160712 ) * liVS[12] +
                   li_Real_s(      -0.000000000022041774031854 ) * liVS[13] +
                   li_Real_s(       0.000000000006957039182342 ) * liVS[14] +
                   li_Real_s(      -0.000000000000940648056671 ) * liVS[15] +
                   li_Real_s(      -0.000000000003543819768056 ) * liVS[16] +
                   li_Real_s(       0.000000000023760539645581 ) * liVS[17] +
                   li_Real_s(      -0.000000000068320562446210 ) * liVS[18] +
                   li_Real_s(       0.000000000108761591067244 ) * liVS[19] +
                   li_Real_s(      -0.000000000103434231452680 ) * liVS[20] +
                   li_Real_s(       0.000000000058826960715695 ) * liVS[21] +
                   li_Real_s(      -0.000000000018558314297224 ) * liVS[22] +
                   li_Real_s(       0.000000000002507836535650 ) * liVS[23] +
                   li_Real_s(       0.000000000005215450643980 ) * liVS[24] +
                   li_Real_s(      -0.000000000035107251942901 ) * liVS[25] +
                   li_Real_s(       0.000000000101178822531980 ) * liVS[26] +
                   li_Real_s(      -0.000000000161243075728507 ) * liVS[27] +
                   li_Real_s(       0.000000000153384041890673 ) * liVS[28] +
                   li_Real_s(      -0.000000000087214745087343 ) * liVS[29] +
                   li_Real_s(       0.000000000027500994924353 ) * liVS[30] +
                   li_Real_s(      -0.000000000003714237232236 ) * liVS[31] +
                   li_Real_s(      -0.000000000004597154314120 ) * liVS[32] +
                   li_Real_s(       0.000000000031075929022651 ) * liVS[33] +
                   li_Real_s(      -0.000000000089790579894578 ) * liVS[34] +
                   li_Real_s(       0.000000000143282250383021 ) * liVS[35] +
                   li_Real_s(      -0.000000000136359492241982 ) * liVS[36] +
                   li_Real_s(       0.000000000077527949128676 ) * liVS[37] +
                   li_Real_s(      -0.000000000024437904299771 ) * liVS[38] +
                   li_Real_s(       0.000000000003299002216102 ) * liVS[39] +
                   li_Real_s(       0.000000000002428549298629 ) * liVS[40] +
                   li_Real_s(      -0.000000000016489972460895 ) * liVS[41] +
                   li_Real_s(       0.000000000047781488777545 ) * liVS[42] +
                   li_Real_s(      -0.000000000076366406666434 ) * liVS[43] +
                   li_Real_s(       0.000000000072725190045711 ) * liVS[44] +
                   li_Real_s(      -0.000000000041352431637268 ) * liVS[45] +
                   li_Real_s(       0.000000000013032273925127 ) * liVS[46] +
                   li_Real_s(      -0.000000000001758691282415 ) * liVS[47] +
                   li_Real_s(      -0.000000000000713364647677 ) * liVS[48] +
                   li_Real_s(       0.000000000004866060141786 ) * liVS[49] +
                   li_Real_s(      -0.000000000014142573169932 ) * liVS[50] +
                   li_Real_s(       0.000000000022643026916432 ) * liVS[51] +
                   li_Real_s(      -0.000000000021581620783159 ) * liVS[52] +
                   li_Real_s(       0.000000000012274759540064 ) * liVS[53] +
                   li_Real_s(      -0.000000000003868175757265 ) * liVS[54] +
                   li_Real_s(       0.000000000000521887759752 ) * liVS[55] +
                   li_Real_s(       0.000000000000090162662007 ) * liVS[56] +
                   li_Real_s(      -0.000000000000617807245876 ) * liVS[57] +
                   li_Real_s(       0.000000000001801017614098 ) * liVS[58] +
                   li_Real_s(      -0.000000000002888768587399 ) * liVS[59] +
                   li_Real_s(       0.000000000002755917440058 ) * liVS[60] +
                   li_Real_s(      -0.000000000001568006932151 ) * liVS[61] +
                   li_Real_s(       0.000000000000494144475482 ) * liVS[62] +
                   li_Real_s(      -0.000000000000066659426219 ) * liVS[63];
        liVC[57] = li_Real_s(       0.000514455781956066338090 ) * liVS[ 0] +
                   li_Real_s(      -0.003601190474163608323810 ) * liVS[ 1] +
                   li_Real_s(       0.010803571423894209324601 ) * liVS[ 2] +
                   li_Real_s(      -0.018005952375225203754283 ) * liVS[ 3] +
                   li_Real_s(       0.018005952376913304802120 ) * liVS[ 4] +
                   li_Real_s(      -0.010803571426922245479751 ) * liVS[ 5] +
                   li_Real_s(       0.003601190475825013712075 ) * liVS[ 6] +
                   li_Real_s(      -0.000514455782277420392568 ) * liVS[ 7] +
                   li_Real_s(      -0.001388888887384559511773 ) * liVS[ 8] +
                   li_Real_s(       0.009722222214678768667517 ) * liVS[ 9] +
                   li_Real_s(      -0.029166666652426932737630 ) * liVS[10] +
                   li_Real_s(       0.048611111099191472995074 ) * liVS[11] +
                   li_Real_s(      -0.048611111108374731315518 ) * liVS[12] +
                   li_Real_s(       0.029166666668976031634442 ) * liVS[13] +
                   li_Real_s(      -0.009722222223844599955922 ) * liVS[14] +
                   li_Real_s(       0.001388888889184543284916 ) * liVS[15] +
                   li_Real_s(       0.002083333330981984587993 ) * liVS[16] +
                   li_Real_s(      -0.014583333324950216447480 ) * liVS[17] +
                   li_Real_s(       0.043749999996405400382038 ) * liVS[18] +
                   li_Real_s(      -0.072916666689692488656505 ) * liVS[19] +
                   li_Real_s(       0.072916666711085376118007 ) * liVS[20] +
                   li_Real_s(      -0.043750000035152634969560 ) * liVS[21] +
                   li_Real_s(       0.014583333346630549323675 ) * liVS[22] +
                   li_Real_s(      -0.002083333335307961664551 ) * liVS[23] +
                   li_Real_s(      -0.002314814813298954332410 ) * liVS[24] +
                   li_Real_s(       0.016203703705319941197338 ) * liVS[25] +
                   li_Real_s(      -0.048611111147117164188458 ) * liVS[26] +
                   li_Real_s(       0.081018518618714480128418 ) * liVS[27] +
                   li_Real_s(      -0.081018518646762779922277 ) * liVS[28] +
                   li_Real_s(       0.048611111198200288030513 ) * liVS[29] +
                   li_Real_s(      -0.016203703734223629667754 ) * liVS[30] +
                   li_Real_s(       0.002314814819167725513244 ) * liVS[31] +
                   li_Real_s(       0.001736111111166900589087 ) * liVS[32] +
                   li_Real_s(      -0.012152777789347868342418 ) * liVS[33] +
                   li_Real_s(       0.036458333395370107199440 ) * liVS[34] +
                   li_Real_s(      -0.060763889026132007764591 ) * liVS[35] +
                   li_Real_s(       0.060763889048472470577611 ) * liVS[36] +
                   li_Real_s(      -0.036458333436303613783735 ) * liVS[37] +
                   li_Real_s(       0.012152777812791412323512 ) * liVS[38] +
                   li_Real_s(      -0.001736111116017388655841 ) * liVS[39] +
                   li_Real_s(      -0.000833333333978260315522 ) * liVS[40] +
                   li_Real_s(       0.005833333344010079545006 ) * liVS[41] +
                   li_Real_s(      -0.017500000046542629056345 ) * liVS[42] +
                   li_Real_s(       0.029166666761542572872123 ) * liVS[43] +
                   li_Real_s(      -0.029166666772330988172124 ) * liVS[44] +
                   li_Real_s(       0.017500000066439896917458 ) * liVS[45] +
                   li_Real_s(      -0.005833333355555872512532 ) * liVS[46] +
                   li_Real_s(       0.000833333336415182398919 ) * liVS[47] +
                   li_Real_s(       0.000231481481824015108972 ) * liVS[48] +
                   li_Real_s(      -0.001620370374663048299024 ) * liVS[49] +
                   li_Real_s(       0.004861111128307893380907 ) * liVS[50] +
                   li_Real_s(      -0.008101851885484852300312 ) * liVS[51] +
                   li_Real_s(       0.008101851888412239799386 ) * liVS[52] +
                   li_Real_s(      -0.004861111133745007681362 ) * liVS[53] +
                   li_Real_s(       0.001620370377861987676305 ) * liVS[54] +
                   li_Real_s(      -0.000231481482513222480701 ) * liVS[55] +
                   li_Real_s(      -0.000028344671260879805708 ) * liVS[56] +
                   li_Real_s(       0.000198412699076934494029 ) * liVS[57] +
                   li_Real_s(      -0.000595238097786278744228 ) * liVS[58] +
                   li_Real_s(       0.000992063496930296884191 ) * liVS[59] +
                   li_Real_s(      -0.000992063497275895433969 ) * liVS[60] +
                   li_Real_s(       0.000595238098432765297754 ) * liVS[61] +
                   li_Real_s(      -0.000198412699462615672863 ) * liVS[62] +
                   li_Real_s(       0.000028344671345714722577 ) * liVS[63];
        liVC[58] = li_Real_s(      -0.000516975307263772831945 ) * liVS[ 0] +
                   li_Real_s(       0.003618827152156112170189 ) * liVS[ 1] +
                   li_Real_s(      -0.010856481460250887172148 ) * liVS[ 2] +
                   li_Real_s(       0.018094135772704578180026 ) * liVS[ 3] +
                   li_Real_s(      -0.018094135777432317213620 ) * liVS[ 4] +
                   li_Real_s(       0.010856481468731540851458 ) * liVS[ 5] +
                   li_Real_s(      -0.003618827156822387348944 ) * liVS[ 6] +
                   li_Real_s(       0.000516975308176724837605 ) * liVS[ 7] +
                   li_Real_s(       0.002212301580405015566377 ) * liVS[ 8] +
                   li_Real_s(      -0.015486111071032253316826 ) * liVS[ 9] +
                   li_Real_s(       0.046458333235732354737912 ) * liVS[10] +
                   li_Real_s(      -0.077430555425377717293500 ) * liVS[11] +
                   li_Real_s(       0.077430555451784233156332 ) * liVS[12] +
                   li_Real_s(      -0.046458333283249851619612 ) * liVS[13] +
                   li_Real_s(       0.015486111097351040744030 ) * liVS[14] +
                   li_Real_s(      -0.002212301585612808096926 ) * liVS[15] +
                   li_Real_s(      -0.004360119032990689147411 ) * liVS[16] +
                   li_Real_s(       0.030520833252879461111595 ) * liVS[17] +
                   li_Real_s(      -0.091562499816971343324923 ) * liVS[18] +
                   li_Real_s(       0.152604166442681071780285 ) * liVS[19] +
                   li_Real_s(      -0.152604166506309202322456 ) * liVS[20] +
                   li_Real_s(       0.091562499931818738807543 ) * liVS[21] +
                   li_Real_s(      -0.030520833316908718702010 ) * liVS[22] +
                   li_Real_s(       0.004360119045800664450141 ) * liVS[23] +
                   li_Real_s(       0.005230379171417304240421 ) * liVS[24] +
                   li_Real_s(      -0.036612654232780414420567 ) * liVS[25] +
                   li_Real_s(       0.109837962782740342082555 ) * liVS[26] +
                   li_Real_s(      -0.183063271418398593182175 ) * liVS[27] +
                   li_Real_s(       0.183063271504780245013677 ) * liVS[28] +
                   li_Real_s(      -0.109837962939127536943218 ) * liVS[29] +
                   li_Real_s(       0.036612654320542586949827 ) * liVS[30] +
                   li_Real_s(      -0.005230379189173914658562 ) * liVS[31] +
                   li_Real_s(      -0.004067460305174880330625 ) * liVS[32] +
                   li_Real_s(       0.028472222165891864636977 ) * liVS[33] +
                   li_Real_s(      -0.085416666571544147545936 ) * liVS[34] +
                   li_Real_s(       0.142361111049476934953262 ) * liVS[35] +
                   li_Real_s(      -0.142361111120639843852231 ) * liVS[36] +
                   li_Real_s(       0.085416666700766794018129 ) * liVS[37] +
                   li_Real_s(      -0.028472222238891595896826 ) * liVS[38] +
                   li_Real_s(       0.004067460320114887895038 ) * liVS[39] +
                   li_Real_s(       0.001994047613835175924990 ) * liVS[40] +
                   li_Real_s(      -0.013958333312982049850826 ) * liVS[41] +
                   li_Real_s(       0.041874999977996658906410 ) * liVS[42] +
                   li_Real_s(      -0.069791666679903704983978 ) * liVS[43] +
                   li_Real_s(       0.069791666715391637021426 ) * liVS[44] +
                   li_Real_s(      -0.041875000042635709962546 ) * liVS[45] +
                   li_Real_s(       0.013958333349745860596514 ) * liVS[46] +
                   li_Real_s(      -0.001994047621447897142288 ) * liVS[47] +
                   li_Real_s(      -0.000561618164562821542418 ) * liVS[48] +
                   li_Real_s(       0.003931327156837742531348 ) * liVS[49] +
                   li_Real_s(      -0.011793981482064383459196 ) * liVS[50] +
                   li_Real_s(       0.019656635817848600034097 ) * liVS[51] +
                   li_Real_s(      -0.019656635827776047786841 ) * liVS[52] +
                   li_Real_s(       0.011793981500202926210319 ) * liVS[53] +
                   li_Real_s(      -0.003931327167225467933420 ) * liVS[54] +
                   li_Real_s(       0.000561618166739462354453 ) * liVS[55] +
                   li_Real_s(       0.000069444444320873599530 ) * liVS[56] +
                   li_Real_s(      -0.000486111110885636184958 ) * liVS[57] +
                   li_Real_s(       0.001458333334132979122733 ) * liVS[58] +
                   li_Real_s(      -0.002430555558689727335686 ) * liVS[59] +
                   li_Real_s(       0.002430555559895151984673 ) * liVS[60] +
                   li_Real_s(      -0.001458333336342092006674 ) * liVS[61] +
                   li_Real_s(       0.000486111112159258451637 ) * liVS[62] +
                   li_Real_s(      -0.000069444444590727183453 ) * liVS[63];
        liVC[59] = li_Real_s(       0.000266479275962155925939 ) * liVS[ 0] +
                   li_Real_s(      -0.001865354932727041059604 ) * liVS[ 1] +
                   li_Real_s(       0.005596064801007345074213 ) * liVS[ 2] +
                   li_Real_s(      -0.009326774672502063934232 ) * liVS[ 3] +
                   li_Real_s(       0.009326774675949521531404 ) * liVS[ 4] +
                   li_Real_s(      -0.005596064807189563006240 ) * liVS[ 5] +
                   li_Real_s(       0.001865354936127624693731 ) * liVS[ 6] +
                   li_Real_s(      -0.000266479276627956673806 ) * liVS[ 7] +
                   li_Real_s(      -0.001406525568782006352375 ) * liVS[ 8] +
                   li_Real_s(       0.009845678987699632220276 ) * liVS[ 9] +
                   li_Real_s(      -0.029537036980005308695052 ) * liVS[10] +
                   li_Real_s(       0.049228394990574839140685 ) * liVS[11] +
                   li_Real_s(      -0.049228395009597726572625 ) * liVS[12] +
                   li_Real_s(       0.029537037014239601179222 ) * liVS[13] +
                   li_Real_s(      -0.009845679006673775657266 ) * liVS[14] +
                   li_Real_s(       0.001406525572544753410753 ) * liVS[15] +
                   li_Real_s(       0.003248181208278522191790 ) * liVS[16] +
                   li_Real_s(      -0.022737268474684338426250 ) * liVS[17] +
                   li_Real_s(       0.068211805467602776298364 ) * liVS[18] +
                   li_Real_s(      -0.113686342505324292084801 ) * liVS[19] +
                   li_Real_s(       0.113686342550596591882694 ) * liVS[20] +
                   li_Real_s(      -0.068211805549360821077087 ) * liVS[21] +
                   li_Real_s(       0.022737268520343242167447 ) * liVS[22] +
                   li_Real_s(      -0.003248181217451705238286 ) * liVS[23] +
                   li_Real_s(      -0.004287918862053485824504 ) * liVS[24] +
                   li_Real_s(       0.030015432059488270088998 ) * liVS[25] +
                   li_Real_s(      -0.090046296241404172322476 ) * liVS[26] +
                   li_Real_s(       0.150077160485185823191046 ) * liVS[27] +
                   li_Real_s(      -0.150077160545851795436789 ) * liVS[28] +
                   li_Real_s(       0.090046296351343452224469 ) * liVS[29] +
                   li_Real_s(      -0.030015432121358882888718 ) * liVS[30] +
                   li_Real_s(       0.004287918874650835203421 ) * liVS[31] +
                   li_Real_s(       0.003506668865650208744000 ) * liVS[32] +
                   li_Real_s(      -0.024546682082246591705132 ) * liVS[33] +
                   li_Real_s(       0.073640046301709422005288 ) * liVS[34] +
                   li_Real_s(      -0.122733410572928192117104 ) * liVS[35] +
                   li_Real_s(       0.122733410622186123184463 ) * liVS[36] +
                   li_Real_s(      -0.073640046391290653327033 ) * liVS[37] +
                   li_Real_s(       0.024546682133060600955465 ) * liVS[38] +
                   li_Real_s(      -0.003506668876140942026076 ) * liVS[39] +
                   li_Real_s(      -0.001772486770630049557340 ) * liVS[40] +
                   li_Real_s(       0.012407407406757097259598 ) * liVS[41] +
                   li_Real_s(      -0.037222222249252014114607 ) * liVS[42] +
                   li_Real_s(       0.062037037117926401019741 ) * liVS[43] +
                   li_Real_s(      -0.062037037142097177511157 ) * liVS[44] +
                   li_Real_s(       0.037222222293372256296529 ) * liVS[45] +
                   li_Real_s(      -0.012407407431993176111229 ) * liVS[46] +
                   li_Real_s(       0.001772486775916681800425 ) * liVS[47] +
                   li_Real_s(       0.000509534832177971352962 ) * liVS[48] +
                   li_Real_s(      -0.003566743828995045176300 ) * liVS[49] +
                   li_Real_s(       0.010700231495535519887774 ) * liVS[50] +
                   li_Real_s(      -0.017833719169499898105080 ) * liVS[51] +
                   li_Real_s(       0.017833719176147816432021 ) * liVS[52] +
                   li_Real_s(      -0.010700231507717587742246 ) * liVS[53] +
                   li_Real_s(       0.003566743836024247049643 ) * liVS[54] +
                   li_Real_s(      -0.000509534833673025433498 ) * liVS[55] +
                   li_Real_s(      -0.000063932980592595889391 ) * liVS[56] +
                   li_Real_s(       0.000447530864639231543145 ) * liVS[57] +
                   li_Real_s(      -0.001342592595007886802083 ) * liVS[58] +
                   li_Real_s(       0.002237654326289389983273 ) * liVS[59] +
                   li_Real_s(      -0.002237654327083241079244 ) * liVS[60] +
                   li_Real_s(       0.001342592596468236872120 ) * liVS[61] +
                   li_Real_s(      -0.000447530865489284210290 ) * liVS[62] +
                   li_Real_s(       0.000063932980776178205407 ) * liVS[63];
        liVC[60] = li_Real_s(      -0.000077160493627145831930 ) * liVS[ 0] +
                   li_Real_s(       0.000540123455700362853360 ) * liVS[ 1] +
                   li_Real_s(      -0.001620370367983209325402 ) * liVS[ 2] +
                   li_Real_s(       0.002700617281247826018120 ) * liVS[ 3] +
                   li_Real_s(      -0.002700617282269033442299 ) * liVS[ 4] +
                   li_Real_s(       0.001620370369811458682863 ) * liVS[ 5] +
                   li_Real_s(      -0.000540123456702315497477 ) * liVS[ 6] +
                   li_Real_s(       0.000077160493822031606115 ) * liVS[ 7] +
                   li_Real_s(       0.000458829364402278325197 ) * liVS[ 8] +
                   li_Real_s(      -0.003211805552786426744272 ) * liVS[ 9] +
                   li_Real_s(       0.009635416663611342039109 ) * liVS[10] +
                   li_Real_s(      -0.016059027779875403862953 ) * liVS[11] +
                   li_Real_s(       0.016059027785278255640034 ) * liVS[12] +
                   li_Real_s(      -0.009635416673331740136654 ) * liVS[13] +
                   li_Real_s(       0.003211805558170041305366 ) * liVS[14] +
                   li_Real_s(      -0.000458829365468343963741 ) * liVS[15] +
                   li_Real_s(      -0.001173941798407223541245 ) * liVS[16] +
                   li_Real_s(       0.008217592594181687484789 ) * liVS[17] +
                   li_Real_s(      -0.024652777795982211417325 ) * liVS[18] +
                   li_Real_s(       0.041087963010750003234328 ) * liVS[19] +
                   li_Real_s(      -0.041087963022983051142312 ) * liVS[20] +
                   li_Real_s(       0.024652777818106371560125 ) * liVS[21] +
                   li_Real_s(      -0.008217592606575063740593 ) * liVS[22] +
                   li_Real_s(       0.001173941800909489296956 ) * liVS[23] +
                   li_Real_s(       0.001679618607418269238263 ) * liVS[24] +
                   li_Real_s(      -0.011757330259958377605134 ) * liVS[25] +
                   li_Real_s(       0.035271990799134990413055 ) * liVS[26] +
                   li_Real_s(      -0.058786651355669426732842 ) * liVS[27] +
                   li_Real_s(       0.058786651371154255429907 ) * liVS[28] +
                   li_Real_s(      -0.035271990827301105686509 ) * liVS[29] +
                   li_Real_s(       0.011757330275936377778079 ) * liVS[30] +
                   li_Real_s(      -0.001679618610714930793115 ) * liVS[31] +
                   li_Real_s(      -0.001455026456764624276374 ) * liVS[32] +
                   li_Real_s(       0.010185185204619395457626 ) * liVS[33] +
                   li_Real_s(      -0.030555555630512844927926 ) * liVS[34] +
                   li_Real_s(       0.050925926070433197301313 ) * liVS[35] +
                   li_Real_s(      -0.050925926082193082033989 ) * liVS[36] +
                   li_Real_s(       0.030555555652044562309211 ) * liVS[37] +
                   li_Real_s(      -0.010185185217012325889496 ) * liVS[38] +
                   li_Real_s(       0.001455026459385722059636 ) * liVS[39] +
                   li_Real_s(       0.000764715609891893077776 ) * liVS[40] +
                   li_Real_s(      -0.005353009273197357684171 ) * liVS[41] +
                   li_Real_s(       0.016059027828278615857016 ) * liVS[42] +
                   li_Real_s(      -0.026765046390189631275724 ) * liVS[43] +
                   li_Real_s(       0.026765046395521244898541 ) * liVS[44] +
                   li_Real_s(      -0.016059027838118741099427 ) * liVS[45] +
                   li_Real_s(       0.005353009278960714489859 ) * liVS[46] +
                   li_Real_s(      -0.000764715611146730349194 ) * liVS[47] +
                   li_Real_s(      -0.000225970018186243126745 ) * liVS[48] +
                   li_Real_s(       0.001581790128504718170172 ) * liVS[49] +
                   li_Real_s(      -0.004745370388051318999434 ) * liVS[50] +
                   li_Real_s(       0.007908950649459373094530 ) * liVS[51] +
                   li_Real_s(      -0.007908950650798274306652 ) * liVS[52] +
                   li_Real_s(       0.004745370390547030969852 ) * liVS[53] +
                   li_Real_s(      -0.001581790129997848876053 ) * liVS[54] +
                   li_Real_s(       0.000225970018522557436480 ) * liVS[55] +
                   li_Real_s(       0.000028935185268523511137 ) * liVS[56] +
                   li_Real_s(      -0.000202546297036853509971 ) * liVS[57] +
                   li_Real_s(       0.000607638891431212455063 ) * liVS[58] +
                   li_Real_s(      -0.001012731486045703038767 ) * liVS[59] +
                   li_Real_s(       0.001012731486190342547471 ) * liVS[60] +
                   li_Real_s(      -0.000607638891704062665371 ) * liVS[61] +
                   li_Real_s(       0.000202546297204210523635 ) * liVS[62] +
                   li_Real_s(      -0.000028935185307676219990 ) * liVS[63];
        liVC[61] = li_Real_s(       0.000012676366845706249498 ) * liVS[ 0] +
                   li_Real_s(      -0.000088734567960707146084 ) * liVS[ 1] +
                   li_Real_s(       0.000266203704000085236703 ) * liVS[ 2] +
                   li_Real_s(      -0.000443672840165559509429 ) * liVS[ 3] +
                   li_Real_s(       0.000443672840286755965078 ) * liVS[ 4] +
                   li_Real_s(      -0.000266203704215788559961 ) * liVS[ 5] +
                   li_Real_s(       0.000088734568077290863590 ) * liVS[ 6] +
                   li_Real_s(      -0.000012676366867762824814 ) * liVS[ 7] +
                   li_Real_s(      -0.000081294091849156736862 ) * liVS[ 8] +
                   li_Real_s(       0.000569058643209471794362 ) * liVS[ 9] +
                   li_Real_s(      -0.001707175930317869193970 ) * liVS[10] +
                   li_Real_s(       0.002845293218073341109942 ) * liVS[11] +
                   li_Real_s(      -0.002845293218629391975016 ) * liVS[12] +
                   li_Real_s(       0.001707175931315745237918 ) * liVS[13] +
                   li_Real_s(      -0.000569058643758485428396 ) * liVS[14] +
                   li_Real_s(       0.000081294091956346167804 ) * liVS[15] +
                   li_Real_s(       0.000223214286360827636468 ) * liVS[16] +
                   li_Real_s(      -0.001562500005251417226418 ) * liVS[17] +
                   li_Real_s(       0.004687500017474205069035 ) * liVS[18] +
                   li_Real_s(      -0.007812500031095116065361 ) * liVS[19] +
                   li_Real_s(       0.007812500032112182704602 ) * liVS[20] +
                   li_Real_s(      -0.004687500019321878992618 ) * liVS[21] +
                   li_Real_s(       0.001562500006294904354720 ) * liVS[22] +
                   li_Real_s(      -0.000223214286573710299355 ) * liVS[23] +
                   li_Real_s(      -0.000340332893712507258965 ) * liVS[24] +
                   li_Real_s(       0.002382330257084740426227 ) * liVS[25] +
                   li_Real_s(      -0.007146990773647027284099 ) * liVS[26] +
                   li_Real_s(       0.011911651291863660506731 ) * liVS[27] +
                   li_Real_s(      -0.011911651292780326555354 ) * liVS[28] +
                   li_Real_s(       0.007146990775351262127624 ) * liVS[29] +
                   li_Real_s(      -0.002382330258094462679952 ) * liVS[30] +
                   li_Real_s(       0.000340332893934677739763 ) * liVS[31] +
                   li_Real_s(       0.000311397708638812276138 ) * liVS[32] +
                   li_Real_s(      -0.002179783961462904758738 ) * liVS[33] +
                   li_Real_s(       0.006539351886386038023580 ) * liVS[34] +
                   li_Real_s(      -0.010898919812445015745106 ) * liVS[35] +
                   li_Real_s(       0.010898919812790138980652 ) * liVS[36] +
                   li_Real_s(      -0.006539351887076839606183 ) * liVS[37] +
                   li_Real_s(       0.002179783961931288337188 ) * liVS[38] +
                   li_Real_s(      -0.000311397708761519675935 ) * liVS[39] +
                   li_Real_s(      -0.000171130953251705242302 ) * liVS[40] +
                   li_Real_s(       0.001197916673299243223508 ) * liVS[41] +
                   li_Real_s(      -0.003593750020899582568040 ) * liVS[42] +
                   li_Real_s(       0.005989583368941778750150 ) * liVS[43] +
                   li_Real_s(      -0.005989583368896064449749 ) * liVS[44] +
                   li_Real_s(       0.003593750020864465693354 ) * liVS[45] +
                   li_Real_s(      -0.001197916673338704929980 ) * liVS[46] +
                   li_Real_s(       0.000171130953280567571495 ) * liVS[47] +
                   li_Real_s(       0.000052358906815307051236 ) * liVS[48] +
                   li_Real_s(      -0.000366512347869878672181 ) * liVS[49] +
                   li_Real_s(       0.001099537043891222511116 ) * liVS[50] +
                   li_Real_s(      -0.001832561740000224204117 ) * liVS[51] +
                   li_Real_s(       0.001832561739925193944334 ) * liVS[52] +
                   li_Real_s(      -0.001099537043772923911034 ) * liVS[53] +
                   li_Real_s(       0.000366512347825504770926 ) * liVS[54] +
                   li_Real_s(      -0.000052358906814201598701 ) * liVS[55] +
                   li_Real_s(      -0.000006889329846380617961 ) * liVS[56] +
                   li_Real_s(       0.000048225308945860532409 ) * liVS[57] +
                   li_Real_s(      -0.000144675926871928092060 ) * liVS[58] +
                   li_Real_s(       0.000241126544804392932420 ) * liVS[59] +
                   li_Real_s(      -0.000241126544787967269506 ) * liVS[60] +
                   li_Real_s(       0.000144675926844819622511 ) * liVS[61] +
                   li_Real_s(      -0.000048225308933972147168 ) * liVS[62] +
                   li_Real_s(       0.000006889329845169347294 ) * liVS[63];
        liVC[62] = li_Real_s(      -0.000001102292774081458163 ) * liVS[ 0] +
                   li_Real_s(       0.000007716049420073507661 ) * liVS[ 1] +
                   li_Real_s(      -0.000023148148265028851197 ) * liVS[ 2] +
                   li_Real_s(       0.000038580247114368220428 ) * liVS[ 3] +
                   li_Real_s(      -0.000038580247116853103387 ) * liVS[ 4] +
                   li_Real_s(       0.000023148148269178472382 ) * liVS[ 5] +
                   li_Real_s(      -0.000007716049421968808255 ) * liVS[ 6] +
                   li_Real_s(       0.000001102292774311403891 ) * liVS[ 7] +
                   li_Real_s(       0.000007440476235672638643 ) * liVS[ 8] +
                   li_Real_s(      -0.000052083333660313526291 ) * liVS[ 9] +
                   li_Real_s(       0.000156250001006275768007 ) * liVS[10] +
                   li_Real_s(      -0.000260416668364911885886 ) * liVS[11] +
                   li_Real_s(       0.000260416668359192231535 ) * liVS[12] +
                   li_Real_s(      -0.000156250000995289574498 ) * liVS[13] +
                   li_Real_s(       0.000052083333653252320830 ) * liVS[14] +
                   li_Real_s(      -0.000007440476233878012997 ) * liVS[15] +
                   li_Real_s(      -0.000021494709146124032877 ) * liVS[16] +
                   li_Real_s(       0.000150462964052693195334 ) * liVS[17] +
                   li_Real_s(      -0.000451388892212083398057 ) * liVS[18] +
                   li_Real_s(       0.000752314820369138310616 ) * liVS[19] +
                   li_Real_s(      -0.000752314820302565261459 ) * liVS[20] +
                   li_Real_s(       0.000451388892093067327187 ) * liVS[21] +
                   li_Real_s(      -0.000150462963987333774984 ) * liVS[22] +
                   li_Real_s(       0.000021494709133207715555 ) * liVS[23] +
                   li_Real_s(       0.000034446649293024770438 ) * liVS[24] +
                   li_Real_s(      -0.000241126545096361096475 ) * liVS[25] +
                   li_Real_s(       0.000723379635349037284318 ) * liVS[26] +
                   li_Real_s(      -0.001205632725549723267486 ) * liVS[27] +
                   li_Real_s(       0.001205632725383285874104 ) * liVS[28] +
                   li_Real_s(      -0.000723379635053685108442 ) * liVS[29] +
                   li_Real_s(       0.000241126544937339808584 ) * liVS[30] +
                   li_Real_s(      -0.000034446649262918068531 ) * liVS[31] +
                   li_Real_s(      -0.000033068783332729874067 ) * liVS[32] +
                   li_Real_s(       0.000231481483369200231404 ) * liVS[33] +
                   li_Real_s(      -0.000694444450142447521242 ) * liVS[34] +
                   li_Real_s(       0.001157407416828561737970 ) * liVS[35] +
                   li_Real_s(      -0.001157407416624434874988 ) * liVS[36] +
                   li_Real_s(       0.000694444449780959663365 ) * liVS[37] +
                   li_Real_s(      -0.000231481483175820441141 ) * liVS[38] +
                   li_Real_s(       0.000033068783296711160014 ) * liVS[39] +
                   li_Real_s(       0.000019014550419887868253 ) * liVS[40] +
                   li_Real_s(      -0.000133101852960357345396 ) * liVS[41] +
                   li_Real_s(       0.000399305558889888171942 ) * liVS[42] +
                   li_Real_s(      -0.000665509264752083429945 ) * liVS[43] +
                   li_Real_s(       0.000665509264614805328732 ) * liVS[44] +
                   li_Real_s(      -0.000399305558646964868941 ) * liVS[45] +
                   li_Real_s(       0.000133101852830779519062 ) * liVS[46] +
                   li_Real_s(      -0.000019014550395955020091 ) * liVS[47] +
                   li_Real_s(      -0.000006062610279276695914 ) * liVS[48] +
                   li_Real_s(       0.000042438271961196641390 ) * liVS[49] +
                   li_Real_s(      -0.000127314815883575876976 ) * liVS[50] +
                   li_Real_s(       0.000212191359780263444530 ) * liVS[51] +
                   li_Real_s(      -0.000212191359731733796348 ) * liVS[52] +
                   li_Real_s(       0.000127314815797735308382 ) * liVS[53] +
                   li_Real_s(      -0.000042438271915493054262 ) * liVS[54] +
                   li_Real_s(       0.000006062610270883995317 ) * liVS[55] +
                   li_Real_s(       0.000000826719583535439652 ) * liVS[56] +
                   li_Real_s(      -0.000005787037085555364336 ) * liVS[57] +
                   li_Real_s(       0.000017361111256372358923 ) * liVS[58] +
                   li_Real_s(      -0.000028935185423261549926 ) * liVS[59] +
                   li_Real_s(       0.000028935185416179351600 ) * liVS[60] +
                   li_Real_s(      -0.000017361111243847994240 ) * liVS[61] +
                   li_Real_s(       0.000005787037078895703314 ) * liVS[62] +
                   li_Real_s(      -0.000000826719582318043242 ) * liVS[63];
        liVC[63] = li_Real_s(       0.000000039367599296174268 ) * liVS[ 0] +
                   li_Real_s(      -0.000000275573195019092779 ) * liVS[ 1] +
                   li_Real_s(       0.000000826719584944077466 ) * liVS[ 2] +
                   li_Real_s(      -0.000001377865974690560124 ) * liVS[ 3] +
                   li_Real_s(       0.000001377865974379526237 ) * liVS[ 4] +
                   li_Real_s(      -0.000000826719584367517385 ) * liVS[ 5] +
                   li_Real_s(       0.000000275573194678167100 ) * liVS[ 6] +
                   li_Real_s(      -0.000000039367599220768007 ) * liVS[ 7] +
                   li_Real_s(      -0.000000275573195340319860 ) * liVS[ 8] +
                   li_Real_s(       0.000001929012367092143788 ) * liVS[ 9] +
                   li_Real_s(      -0.000005787037100368583144 ) * liVS[10] +
                   li_Real_s(       0.000009645061831904932223 ) * liVS[11] +
                   li_Real_s(      -0.000009645061829056682294 ) * liVS[12] +
                   li_Real_s(       0.000005787037095193399030 ) * liVS[13] +
                   li_Real_s(      -0.000001929012364155480559 ) * liVS[14] +
                   li_Real_s(       0.000000275573194730584887 ) * liVS[15] +
                   li_Real_s(       0.000000826719586357191002 ) * liVS[16] +
                   li_Real_s(      -0.000005787037103760761209 ) * liVS[17] +
                   li_Real_s(       0.000017361111308213153764 ) * liVS[18] +
                   li_Real_s(      -0.000028935185506414442321 ) * liVS[19] +
                   li_Real_s(       0.000028935185496548348241 ) * liVS[20] +
                   li_Real_s(      -0.000017361111290434918147 ) * liVS[21] +
                   li_Real_s(       0.000005787037093856737841 ) * liVS[22] +
                   li_Real_s(      -0.000000826719584365301547 ) * liVS[23] +
                   li_Real_s(      -0.000001377865977234837033 ) * liVS[24] +
                   li_Real_s(       0.000009645061839512975379 ) * liVS[25] +
                   li_Real_s(      -0.000028935185512954241405 ) * liVS[26] +
                   li_Real_s(       0.000048225309174876252372 ) * liVS[27] +
                   li_Real_s(      -0.000048225309157152982417 ) * liVS[28] +
                   li_Real_s(       0.000028935185481146050206 ) * liVS[29] +
                   li_Real_s(      -0.000009645061821961261784 ) * liVS[30] +
                   li_Real_s(       0.000001377865973768067129 ) * liVS[31] +
                   li_Real_s(       0.000001377865976853990691 ) * liVS[32] +
                   li_Real_s(      -0.000009645061836877742378 ) * liVS[33] +
                   li_Real_s(       0.000028935185504664333338 ) * liVS[34] +
                   li_Real_s(      -0.000048225309160238662457 ) * liVS[35] +
                   li_Real_s(       0.000048225309141831599745 ) * liVS[36] +
                   li_Real_s(      -0.000028935185471700785811 ) * liVS[37] +
                   li_Real_s(       0.000009645061818788364455 ) * liVS[38] +
                   li_Real_s(      -0.000001377865973321090807 ) * liVS[39] +
                   li_Real_s(      -0.000000826719585760752752 ) * liVS[40] +
                   li_Real_s(       0.000005787037099671887533 ) * liVS[41] +
                   li_Real_s(      -0.000017361111295279231710 ) * liVS[42] +
                   li_Real_s(       0.000028935185483340014617 ) * liVS[43] +
                   li_Real_s(      -0.000028935185472149709885 ) * liVS[44] +
                   li_Real_s(       0.000017361111275266052228 ) * liVS[45] +
                   li_Real_s(      -0.000005787037088728129528 ) * liVS[46] +
                   li_Real_s(       0.000000826719583639865262 ) * liVS[47] +
                   li_Real_s(       0.000000275573195111629434 ) * liVS[48] +
                   li_Real_s(      -0.000001929012365569338214 ) * liVS[49] +
                   li_Real_s(       0.000005787037095444358799 ) * liVS[50] +
                   li_Real_s(      -0.000009645061822791267675 ) * liVS[51] +
                   li_Real_s(       0.000009645061819087351839 ) * liVS[52] +
                   li_Real_s(      -0.000005787037088826192227 ) * liVS[53] +
                   li_Real_s(       0.000001929012361959605242 ) * liVS[54] +
                   li_Real_s(      -0.000000275573194416145928 ) * liVS[55] +
                   li_Real_s(      -0.000000039367599279396240 ) * liVS[56] +
                   li_Real_s(       0.000000275573194926447492 ) * liVS[57] +
                   li_Real_s(      -0.000000826719584600192254 ) * liVS[58] +
                   li_Real_s(       0.000001377865973917845847 ) * liVS[59] +
                   li_Real_s(      -0.000001377865973400771196 ) * liVS[60] +
                   li_Real_s(       0.000000826719583676820356 ) * liVS[61] +
                   li_Real_s(      -0.000000275573194423754825 ) * liVS[62] +
                   li_Real_s(       0.000000039367599182949680 ) * liVS[63];

        /* Prepare interpolated value computation */
        liTX1 = ( liX + li_Real_s( 3.0 ) ) - liPX; 
        liTX2 = liTX1 * liTX1;
        liTX3 = liTX1 * liTX2;
        liTX4 = liTX1 * liTX3;
        liTX5 = liTX1 * liTX4;
        liTX6 = liTX1 * liTX5;
        liTX7 = liTX1 * liTX6;
        liTY1 = ( liY + li_Real_s( 3.0 ) ) - liPY; 
        liTY2 = liTY1 * liTY1; 
        liTY3 = liTY1 * liTY2; 
        liTY4 = liTY1 * liTY3; 
        liTY5 = liTY1 * liTY4; 
        liTY6 = liTY1 * liTY5; 
        liTY7 = liTY1 * liTY6;

        /* Compute interpolated value */
        liIV = liVC[ 0]                 + 
               liVC[ 1] * liTY1         + 
               liVC[ 2] * liTY2         + 
               liVC[ 3] * liTY3         + 
               liVC[ 4] * liTY4         + 
               liVC[ 5] * liTY5         + 
               liVC[ 6] * liTY6         + 
               liVC[ 7] * liTY7         +
               liVC[ 8] * liTX1         + 
               liVC[ 9] * liTY1 * liTX1 + 
               liVC[10] * liTY2 * liTX1 + 
               liVC[11] * liTY3 * liTX1 + 
               liVC[12] * liTY4 * liTX1 + 
               liVC[13] * liTY5 * liTX1 + 
               liVC[14] * liTY6 * liTX1 + 
               liVC[15] * liTY7 * liTX1 +
               liVC[16] * liTX2         + 
               liVC[17] * liTY1 * liTX2 + 
               liVC[18] * liTY2 * liTX2 + 
               liVC[19] * liTY3 * liTX2 + 
               liVC[20] * liTY4 * liTX2 + 
               liVC[21] * liTY5 * liTX2 + 
               liVC[22] * liTY6 * liTX2 + 
               liVC[23] * liTY7 * liTX2 +
               liVC[24] * liTX3         + 
               liVC[25] * liTY1 * liTX3 + 
               liVC[26] * liTY2 * liTX3 + 
               liVC[27] * liTY3 * liTX3 + 
               liVC[28] * liTY4 * liTX3 + 
               liVC[29] * liTY5 * liTX3 + 
               liVC[30] * liTY6 * liTX3 + 
               liVC[31] * liTY7 * liTX3 +
               liVC[32] * liTX4         + 
               liVC[33] * liTY1 * liTX4 + 
               liVC[34] * liTY2 * liTX4 + 
               liVC[35] * liTY3 * liTX4 + 
               liVC[36] * liTY4 * liTX4 + 
               liVC[37] * liTY5 * liTX4 + 
               liVC[38] * liTY6 * liTX4 + 
               liVC[39] * liTY7 * liTX4 +
               liVC[40] * liTX5         + 
               liVC[41] * liTY1 * liTX5 + 
               liVC[42] * liTY2 * liTX5 + 
               liVC[43] * liTY3 * liTX5 + 
               liVC[44] * liTY4 * liTX5 + 
               liVC[45] * liTY5 * liTX5 + 
               liVC[46] * liTY6 * liTX5 + 
               liVC[47] * liTY7 * liTX5 +
               liVC[48] * liTX6         + 
               liVC[49] * liTY1 * liTX6 + 
               liVC[50] * liTY2 * liTX6 + 
               liVC[51] * liTY3 * liTX6 + 
               liVC[52] * liTY4 * liTX6 + 
               liVC[53] * liTY5 * liTX6 + 
               liVC[54] * liTY6 * liTX6 + 
               liVC[55] * liTY7 * liTX6 +  
               liVC[56] * liTX7         + 
               liVC[57] * liTY1 * liTX7 + 
               liVC[58] * liTY2 * liTX7 + 
               liVC[59] * liTY3 * liTX7 + 
               liVC[60] * liTY4 * liTX7 + 
               liVC[61] * liTY5 * liTX7 + 
               liVC[62] * liTY6 * liTX7 + 
               liVC[63] * liTY7 * liTX7;

        /* Verify interpolated value */
        liIV = ( liIV < li_Real_s(   0.0 ) ) ? li_Real_s(   0.0 ) : liIV; 
        liIV = ( liIV > li_Real_s( 255.0 ) ) ? li_Real_s( 255.0 ) : liIV;

        /* Return interpolated value */
        return( liIV );

    }

