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
    //! @author Nils Hamel (nils.hamel@foxel.ch)
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

    //! This function performe an order four fast bicubic interpolation of bitmap 
    //! pixels based on bm bitmap. The value of floating point coordinates have to 
    //! be in the [0,bm_w-1[ × [0,bm_h-1[ range. This last condition is not verified
    //! by the function.
    //! 
    //! @param bm Pointer to bitmap array
    //! @param bm_w Bitmap width
    //! @param bm_h Bitmap height
    //! @param bm_d Bitmap number of chromatic layer
    //! @param bm_c Bitmap interpolated layer
    //! @param bm_x Interpolated point position x
    //! @param bm_y Interpolated point position y    
    //! @return Interpolated value

    inter_C8_t inter_bicubicf(

        inter_C8_t *  bm, 
        inter_Index_t bm_w,
        inter_Index_t bm_h,
        inter_Index_t bm_d, 
        inter_Index_t bm_c,
        inter_Real_t  bm_x,
        inter_Real_t  bm_y

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

