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

    //! @file   inter-bicubic.h
    //! @author Nils Hamel (n.hamel@foxel.ch)
    //! 
    //! Bicubic interpolation methods

/* 
    Header - Include guard
 */

    # ifndef __LIBINTER_BICUBIC__
    # define __LIBINTER_BICUBIC__

/* 
    Header - C/C++ compatibility
 */

    # ifdef __cplusplus
    extern "C" {
    # endif

/* 
    Header - Includes
 */

    # include "inter.h"

/* 
    Header - Preprocessor definitions
 */

/* 
    Header - Preprocessor macros
 */

/* 
    Header - Typedefs
 */

/* 
    Header - Structures
 */

/* 
    Header - Function prototypes
 */

    //! Fast bicubic interpolation method

    //! This function performe an order four fast bicubic interpolation of bitmap 
    //! pixels based on bm bitmap. The value of floating point coordinates have to 
    //! be in the [0,liWidth-1[ × [0,liHeight-1[ range. This last condition is not
    //! verified by the function.
    //! 
    //! @param liBytes Pointer to bitmap array
    //! @param liWidth Bitmap width
    //! @param liHeight Bitmap height
    //! @param liLayer Bitmap number of chromatic layer
    //! @param liChannel Bitmap interpolated layer
    //! @param liX Interpolated point position x (floating point)
    //! @param liY Interpolated point position y (floating point)
    //! @return Interpolated value

    inter_C8_t inter_bicubicf(

        inter_C8_t * liBytes, 
        inter_Size_t liWidth,
        inter_Size_t liHeight,
        inter_Size_t liLayer, 
        inter_Size_t liChannel,
        inter_Real_t liX,
        inter_Real_t liY

    );

/* 
    Header - C/C++ compatibility
 */

    # ifdef __cplusplus
    } 
    # endif

/*
    Header - Include guard
 */

    # endif

