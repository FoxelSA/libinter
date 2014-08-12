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

    //! @file   inter.h
    //! @author Nils Hamel (n.hamel@foxel.ch)
    //! 
    //! Library common header file

/* 
    Header - Include guard
 */

    # ifndef __LIBINTER_INTER__
    # define __LIBINTER_INTER__

/* 
    Header - C/C++ compatibility
 */

    # ifdef __cplusplus
    extern "C" {
    # endif

/* 
    Header - Includes
 */

    # include <math.h>
    # include <stdint.h>

/* 
    Header - Preprocessor definitions
 */

/* 
    Header - Preprocessor macros
 */

/* 
    Header - Typedefs
 */

    /* Define pixel component type */
    typedef uint8_t inter_C8_t;

    /* Define general index */
    typedef int64_t inter_Size_t;

    /* Define general enumeration */
    typedef int64_t inter_Enum_t;

    /* Define floating type */
    typedef double inter_Real_t;

    /* General interpolation method prototype */
    typedef inter_C8_t ( * inter_Method_t ) ( inter_C8_t * , inter_Size_t , inter_Size_t , inter_Size_t , inter_Size_t , inter_Real_t , inter_Real_t );

    /* Define literal suffix */
    # define inter_C8_s  ( x )  UINT8_C( x )
    # define inter_Size_s( x )  INT64_C( x )
    # define inter_Enum_s( x )  INT64_C( x )
    # define inter_Real_s( x )  ( x )

    /* Define formated output specifiers */
    # define inter_C8_p         PRIu8
    # define inter_Size_p       PRId64
    # define inter_Enum_p       PRId64
    # define inter_Real_p       "lf"

    /* Define formated input specifiers */
    # define inter_C8_i         SCNu8
    # define inter_Size_i       SCNu64
    # define inter_Enum_i       SCNu64
    # define inter_Real_i       "lf"

/* 
    Header - Structures
 */

/* 
    Header - Function prototypes
 */

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

