/*
   libit - Library for basic source and channel coding functions
   Copyright (C) 2005-2005 Vivien Chappelier, Herve Jegou

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public
   License as published by the Free Software Foundation; either
   version 2 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU Library General Public
   License along with this library; if not, write to the Free
   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/ 

/*
  Random numbers generators

  This Libit code is: 
  Copyright (C) 2005-2006 Vivien Chappelier, Fran�ois Cayre  
  
  Based on the mt19937 implementation (LGPL)
    Copyright (C) 1997, 1999 Makoto Matsumoto and Takuji Nishimura.
    http://www.math.keio.ac.jp/matumoto/emt.html
  REFERENCE
    M. Matsumoto and T. Nishimura,
    "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform
     Pseudo-Random Number Generator",                             
    ACM Transactions on Modeling and Computer Simulation,
    Vol. 8, No. 1, January 1998, pp 3--30.

  Uses ZIGGURAT method for GAUSSIAN pseudo-random numbers generation 
     Copyright (C) 2006 Fran�ois Cayre

  Our Gaussian PRNG uses 256 ziggurat levels and tries to stick 
  as much as possible to the original paper. In particular, 
  compared to the GSL implementation of J. Voss, we do not model 
  the tail of the Gaussian with an exponential distribution. 
  Informal tests on an x86_64 machine do not reveal significant 
  speed improvement of one implementation over the other. 

  REFERENCE 
     George Marsagli, Wai Wan Tsang, 
     "The Ziggurat Method for Generating Random Variables",
     Journal of Statistical Software, vol. 5 (2000), no. 8.


    TODO: 
    - Add exponential PRNG based on Ziggurat

    Changelog (06/16/06) 
    - Replaced Box-Muller based it_randn() function by Ziggurat method 

    Changelog (06/15/06) 
    - Now based on: 
    http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/MTARCOK/mt19937ar-cok.c 
    - Provides more functions for various real output 

*/


#include <math.h>
#ifdef WIN32
#include <windows.h> /* sigh... */
#else
#include <sys/time.h>
#endif
#include <stdlib.h>
#include "../include/vec.h"
#include "../include/random.h"

typedef unsigned int it_uint32_t; /* Solaris and Windows don't have stdint... */

/* MT19937 related stuff */

/* MT199737 constants */  
#define MSB 0x80000000UL /* most significant bit */
#define LSB 0x7fffffffUL /* least significant bits */
#define N          624
#define M          397
#define A   0x9908b0dfUL  
#define MIXBITS(u,v) ( ((u) & MSB) | ((v) & LSB) )
#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? A : 0UL))


/* some random default state (generated from a random seed) */
#define MT199737_STATE_INITIALIZER { 0x7036c3e9UL, 0x67b6b695UL, 0x0b2ac651UL, 0xe3a3ebddUL, 0x891b18f9UL, 0x43f0a465UL, 0x4471c5e1UL, 0x5c08e22dUL, 0x08c66709UL, 0x8124f735UL, 0x425ca671UL, 0x7f16057dUL, 0xe866be19UL, 0x10067f05UL, 0xa8abf801UL, 0xc54ea5cdUL, 0xa4332e29UL, 0x447d0bd5UL, 0x3d914a91UL, 0xf8b3131dUL, 0x527bc739UL, 0x1d756da5UL, 0xfb3f2e21UL, 0x95109d6dUL, 0xd3b99949UL, 0x91b17475UL, 0x227932b1UL, 0xf55194bdUL, 0x2b9eb459UL, 0xb197f045UL, 0xde23e841UL, 0x9dcd490dUL, 0x13262869UL, 0x6e04b115UL, 0x07d4ded1UL, 0xf3980a5dUL, 0xd3a40579UL, 0xd41886e5UL, 0x9d62a661UL, 0xb0d328adUL, 0x7ad55b89UL, 0x8e0941b5UL, 0x7774cef1UL, 0x65fcf3fdUL, 0x77f03a99UL, 0x78f1b185UL, 0xd113e881UL, 0x5840bc4dUL, 0xb1b3b2a9UL, 0x1fa1a655UL, 0x30398311UL, 0x0cc6d19dUL, 0x2577d3b9UL, 0xea6df025UL, 0x3f602ea1UL, 0xd10483edUL, 0x1f3dadc9UL, 0xd4005ef5UL, 0x28137b31UL, 0x900c233dUL, 0x1abf50d9UL, 0x7327c2c5UL, 0xfe7ff8c1UL, 0x44dcff8dUL, 0x5d7fcce9UL, 0x29a7eb95UL, 0xce033751UL, 0x59b368ddUL, 0x59db31f9UL, 0x4809a965UL, 0xd6bbc6e1UL, 0x4458af2dUL, 0xeb169009UL, 0xf66acc35UL, 0x8c193771UL, 0xaf73227dUL, 0x7a6ff719UL, 0xb24e2405UL, 0xe46c1901UL, 0x80d612cdUL, 0x1d2e7729UL, 0x316b80d5UL, 0x8975fb91UL, 0x0cd1d01dUL, 0xbbb22039UL, 0x797fb2a5UL, 0x79f96f21UL, 0xc683aa6dUL, 0x51840249UL, 0xfd1c8975UL, 0xac4a03b1UL, 0xbd25f1bdUL, 0x56662d59UL, 0x8d78d545UL, 0x41dc4941UL, 0x365ff60dUL, 0x6063b169UL, 0xf1406615UL, 0xdbd5cfd1UL, 0xb596075dUL, 0x0ee09e79UL, 0xf0640be5UL, 0xa09d2761UL, 0xc03975adUL, 0x4eaa0489UL, 0xa4e996b5UL, 0x8269dff1UL, 0xaf1890fdUL, 0x0705f399UL, 0xe0bbd685UL, 0x56d48981UL, 0xdcaea94dUL, 0x3fc37ba9UL, 0x787a9b55UL, 0x4f66b411UL, 0x80740e9dUL, 0xd04aacb9UL, 0x434ab525UL, 0x632aefa1UL, 0x872e10edUL, 0xa7ac96c9UL, 0x9fa5f3f5UL, 0x393ccc31UL, 0xb83f003dUL, 0xbdb349d9UL, 0x4d2b27c5UL, 0x2458d9c1UL, 0x77f62c8dUL, 0xbcf1d5e9UL, 0x6b6e2095UL, 0xbf6ca851UL, 0x76dfe5ddUL, 0x75d44af9UL, 0x6dc7ae65UL, 0xbb26c7e1UL, 0x9e157c2dUL, 0x2aafb909UL, 0xd425a135UL, 0x6c86c871UL, 0x888d3f7dUL, 0xc4d23019UL, 0x78dac905UL, 0xac6d3a01UL, 0xd96a7fcdUL, 0x0292c029UL, 0x436ef5d5UL, 0x982bac91UL, 0xbf4d8d1dUL, 0xae617939UL, 0x106ef7a5UL, 0xc314b021UL, 0xf4a3b76dUL, 0xeed76b49UL, 0x9e3c9e75UL, 0x690bd4b1UL, 0x8cf74ebdUL, 0xbfc6a659UL, 0x4edeba45UL, 0x3215aa41UL, 0xdf3fa30dUL, 0xa44a3a69UL, 0x8ed11b15UL, 0x16e7c0d1UL, 0xdd31045dUL, 0xa1d63779UL, 0xb0d490e5UL, 0xf678a861UL, 0x278cc2adUL, 0x9447ad89UL, 0x0ebeebb5UL, 0x6c8ff0f1UL, 0x2f712dfdUL, 0xeaf4ac99UL, 0x3f4afb85UL, 0x79562a81UL, 0xb4a9964dUL, 0xdebc44a9UL, 0x30e89055UL, 0x89e4e511UL, 0xf0fe4b9dUL, 0x311685b9UL, 0xf98c7a25UL, 0x71d6b0a1UL, 0xc0849dedUL, 0x84247fc9UL, 0x2b8088f5UL, 0xe5d71d31UL, 0x16eedd3dUL, 0x5bc042d9UL, 0x7f338cc5UL, 0x0732bac1UL, 0x11dc598dUL, 0xd78cdee9UL, 0xa2095595UL, 0x90671951UL, 0xf82962ddUL, 0x360663f9UL, 0xfa2ab365UL, 0x32b2c8e1UL, 0x763f492dUL, 0x3091e209UL, 0x2f557635UL, 0xb4a55971UL, 0x67645c7dUL, 0x408d6919UL, 0x48ac6e05UL, 0x61af5b01UL, 0x7c0beccdUL, 0xdd600929UL, 0x2f876ad5UL, 0x5ab25d91UL, 0x0d264a1dUL, 0xc389d239UL, 0x67433ca5UL, 0x5790f121UL, 0x6c70c46dUL, 0x54b3d449UL, 0xca11b375UL, 0x69bea5b1UL, 0x01c5abbdUL, 0x20c01f59UL, 0x1ac99f45UL, 0x4fd00b41UL, 0x856c500dUL, 0xa7d9c369UL, 0x3bb6d015UL, 0xea0ab1d1UL, 0xa769015dUL, 0x6584d079UL, 0xda6a15e5UL, 0x5ff52961UL, 0x73cd0fadUL, 0x34ae5689UL, 0x608940b5UL, 0x86e701f1UL, 0xc406cafdUL, 0x1cbc6599UL, 0xf99f2085UL, 0x1998cb81UL, 0x0d31834dUL, 0x979e0da9UL, 0x7deb8555UL, 0x50b41611UL, 0xdb65889dUL, 0x60db5eb9UL, 0x12333f25UL, 0x6c6371a1UL, 0x4a082aedUL, 0xdda568c9UL, 0x4c901df5UL, 0xbee26e31UL, 0xc91bba3dUL, 0x2de63bd9UL, 0xae40f1c5UL, 0xc80d9bc1UL, 0x7f8f868dUL, 0xf650e7e9UL, 0x42798a95UL, 0xf1f28a51UL, 0x9a8fdfddUL, 0xf3717cf9UL, 0x3232b865UL, 0x7e5fc9e1UL, 0xd9d6162dUL, 0x65bd0b09UL, 0x1cfa4b35UL, 0x3574ea71UL, 0xa8f8797dUL, 0x66a1a219UL, 0x06c31305UL, 0x65327c01UL, 0x15ba59cdUL, 0x36965229UL, 0xaab4dfd5UL, 0xc20a0e91UL, 0xf35c071dUL, 0x942b2b39UL, 0x02fc81a5UL, 0xb86e3221UL, 0x7aead16dUL, 0x2c193d49UL, 0xd59bc875UL, 0xbf6276b1UL, 0xb89108bdUL, 0x32529859UL, 0x16398445UL, 0x3c0b6c41UL, 0x15e5fd0dUL, 0x34124c69UL, 0xecf18515UL, 0x863ea2d1UL, 0x513dfe5dUL, 0x32ec6979UL, 0x32249ae5UL, 0x9e12aa61UL, 0x31fa5cadUL, 0x18ddff89UL, 0x2f4895b5UL, 0x226f12f1UL, 0x49d967fdUL, 0x955d1e99UL, 0x74b84585UL, 0x189c6c81UL, 0x1346704dUL, 0x7368d6a9UL, 0x94837a55UL, 0x14d44711UL, 0xbca9c59dUL, 0x789937b9UL, 0x923f0425UL, 0x53d132a1UL, 0xf0b8b7edUL, 0xdd2f51c9UL, 0xd7d4b2f5UL, 0x555ebf31UL, 0xebc5973dUL, 0x6d2534d9UL, 0x7f5356c5UL, 0x87e97cc1UL, 0x2e0fb38dUL, 0x623df0e9UL, 0xc1bebf95UL, 0x950efb51UL, 0x1b135cddUL, 0x071595f9UL, 0x5adfbd65UL, 0xdf2dcae1UL, 0xd5d9e32dUL, 0x33313409UL, 0xb2142035UL, 0xbff57b71UL, 0xaa49967dUL, 0xb00edb19UL, 0x981eb805UL, 0x17f69d01UL, 0x5375c6cdUL, 0x97359b29UL, 0x69f754d5UL, 0xbf32bf91UL, 0x6eeec41dUL, 0xb9458439UL, 0x689ac6a5UL, 0x66ac7321UL, 0x6d11de6dUL, 0x1e07a649UL, 0x15dadd75UL, 0x7af747b1UL, 0x4e5965bdUL, 0xad7e1159UL, 0x662e6945UL, 0x97c7cd41UL, 0x7dacaa0dUL, 0x11f3d569UL, 0x97813a15UL, 0x1c8393d1UL, 0x17affb5dUL, 0xe30d0279UL, 0x7d041fe5UL, 0x71d12b61UL, 0xef14a9adUL, 0x29d6a889UL, 0x0ffceab5UL, 0x902823f1UL, 0x9de904fdUL, 0x4dd6d799UL, 0x15966a85UL, 0x57610d81UL, 0xf3e85d4dUL, 0x7b1c9fa9UL, 0xa9b06f55UL, 0x47457811UL, 0x11cb029dUL, 0x915010b9UL, 0x7eafc925UL, 0x291ff3a1UL, 0x819644edUL, 0xabc23ac9UL, 0xa24e47f5UL, 0x3a4c1031UL, 0x9bec743dUL, 0x527d2dd9UL, 0x976abbc5UL, 0x67c65dc1UL, 0x8a5ce08dUL, 0x6453f9e9UL, 0x94d8f495UL, 0x2abc6c51UL, 0x36b3d9ddUL, 0xc9f2aef9UL, 0xb931c265UL, 0x961ccbe1UL, 0x774ab02dUL, 0x01ee5d09UL, 0x03a2f535UL, 0x25270c71UL, 0xc857b37dUL, 0x95d51419UL, 0xe1bf5d05UL, 0xdafbbe01UL, 0xe23e33cdUL, 0x883de429UL, 0x224ec9d5UL, 0x432c7091UL, 0x7cde811dUL, 0xcbd8dd39UL, 0x1d1e0ba5UL, 0xe34bb421UL, 0x8fe5eb6dUL, 0xd37f0f49UL, 0xdfcef275UL, 0xad7d18b1UL, 0x601ec2bdUL, 0x4b428a59UL, 0x2fa84e45UL, 0x04052e41UL, 0xa9c0570dUL, 0x0a7e5e69UL, 0x3065ef15UL, 0xddd984d1UL, 0x37bef85dUL, 0x4ee69b79UL, 0x8008a4e5UL, 0x9c30ac61UL, 0x381bf6adUL, 0x50985189UL, 0x97a63fb5UL, 0x211234f1UL, 0x9d35a1fdUL, 0x3f299099UL, 0x41398f85UL, 0xb6e6ae81UL, 0xdc174a4dUL, 0xb7b968a9UL, 0xf2726455UL, 0x5907a911UL, 0x57c93f9dUL, 0xc3ffe9b9UL, 0xdc858e25UL, 0xed4fb4a1UL, 0xc9a0d1edUL, 0x725e23c9UL, 0x80fcdcf5UL, 0xfeaa6131UL, 0xf690513dUL, 0x16ee26d9UL, 0x9b8720c5UL, 0x88a43ec1UL, 0x01770d8dUL, 0x459302e9UL, 0x30c82995UL, 0x63fadd51UL, 0xaa7156ddUL, 0x9508c7f9UL, 0x9228c765UL, 0xe42ccce1UL, 0xcb287d2dUL, 0x3af48609UL, 0x26a6ca35UL, 0x36099d71UL, 0x6022d07dUL, 0x90f44d19UL, 0xc8a50205UL, 0x0f41df01UL, 0x6f13a0cdUL, 0x92af2d29UL, 0x88bb3ed5UL, 0x3ef72191UL, 0x1a2b3e1dUL, 0x64e53639UL, 0xa58650a5UL, 0xaf4bf521UL, 0x3066f86dUL, 0xf57f7849UL, 0x88780775UL, 0x67f3e9b1UL, 0x8ae11fbdUL, 0xc4a00359UL, 0x97a73345UL, 0x21c38f41UL, 0x8721040dUL, 0xe6b1e769UL, 0xac9fa415UL, 0xfb4075d1UL, 0xee6af55dUL, 0x4f793479UL, 0x003229e5UL, 0xde312d61UL, 0x9a1043adUL, 0x7622fa89UL, 0x5b4494b5UL, 0x262d45f1UL, 0x24bf3efdUL, 0x62554999UL, 0x5ca1b485UL, 0x182d4f81UL, 0xf8d3374dUL, 0x323f31a9UL, 0xa3c95955UL, 0xbb1ada11UL, 0x0ba47c9dUL, 0x29a8c2b9UL, 0xb0c05325UL, 0xa16075a1UL, 0x95d85eedUL, 0x5a030cc9UL, 0x48e071f5UL, 0x3379b231UL, 0x18b12e3dUL, 0xf3781fd9UL, 0x30a885c5UL, 0x0b831fc1UL, 0x005e3a8dUL, 0x4efb0be9UL, 0x0a8c5e95UL, 0xf1ca4e51UL, 0x334bd3ddUL, 0xc157e0f9UL, 0x2ac4cc65UL, 0x0a5dcde1UL, 0xde734a2dUL, 0x4743af09UL, 0x301f9f35UL, 0xc39d2e71UL, 0xceaaed7dUL, 0x1a6c8619UL, 0x31cfa705UL, 0x15c90001UL, 0xa6f60dcdUL, 0x3f897629UL, 0x523cb3d5UL, 0xa392d291UL, 0x43d4fb1dUL, 0x1d6a8f39UL, 0x86d395a5UL, 0x4bad3621UL, 0x9b95056dUL, 0x2d08e149UL, 0x64d61c75UL, 0xbb5bbab1UL, 0x6ba07cbdUL, 0xd2967c59UL, 0xc32b1845UL, 0x9202f041UL, 0x02ceb10dUL, 0x6f8e7069UL, 0x012e5915UL, 0xa5b866d1UL, 0x78b3f25dUL, 0xbdc4cd79UL, 0xc280aee5UL, 0xf8d2ae61UL, 0xa1f190adUL, 0x8376a389UL, 0xefd7e9b5UL, 0xf07956f1UL, 0x1185dbfdUL, 0xb05a0299UL, 0xccced985UL, 0x5c34f081UL, 0x771c244dUL, 0xf3adfaa9UL, 0xf2b54e55UL, 0xde7f0b11UL, 0xaa5cb99dUL, 0xdb4a9bb9UL, 0x00601825UL, 0x465236a1UL, 0xb33cebedUL, 0x8bb0f5c9UL, 0xcef906f5UL, 0x69ba0331UL, 0x1f4f0b3dUL, 0x211b18d9UL, 0xfbceeac5UL, 0x116300c1UL, 0xf412678dUL, 0xc98c14e9UL, 0x97259395UL, 0x852abf51UL, 0x8e4350ddUL, 0xa7dff9f9UL, 0xc805d165UL, 0x49afcee1UL, 0xbe2b172dUL, 0x8fdbd809UL, 0x350d7435UL, 0x9ee1bf71UL, 0x70f00a7dUL, 0xab3dbf19UL, 0x023f4c05UL, 0x4f912101UL, 0x36e57acdUL, 0x17ccbf29UL, 0x33d328d5UL, 0x61ff8391UL, 0xf6dbb81dUL, 0x8e68e839UL, 0x4605daa5UL, 0x396f7721UL, 0x1e70126dUL, 0x231b4a49UL, 0xc9e93175UL, 0xb8b48bb1UL, 0x9f5cd9bdUL, 0x2e25f559UL, 0xd733fd45UL, 0xf5c35141UL, 0x09c95e0dUL, 0x6e13f969UL, 0x23120e15UL, 0x0e4157d1UL, 0x1399ef5dUL, 0x72c96679UL, 0x8bf433e5UL, 0xad152f61UL, 0xdcbfddadUL, 0x61934c89UL, 0xea603eb5UL, 0xd0f667f1UL, 0x408978fdUL, 0x2237bb99UL, 0xf6c0fe85UL, 0x63fd9181UL, 0x83f2114dUL }

static it_uint32_t state[N] = MT199737_STATE_INITIALIZER; /* state vector  */
static int left = 1; 
static int initf = 0; 
static it_uint32_t *next; 

static int idx = 0;                                    /* current index */

/* Ziggurat related stuff */


#define ZIGLEVELS 256

#define ZIGR 3.6554204190269

#define ZIGRINV 0.2735663440503

#define ZIGV 0.0049287609634869

static const unsigned long int zkn[ZIGLEVELS] = {
3995918566, 3230394625, 3660926751, 3844040550, 
3945094684, 4009015195, 4053037110, 4085175722, 
4109658338, 4128923088, 4144473829, 4157288020, 
4168028249, 4177159436, 4185017354, 4191850516, 
4197846767, 4203150717, 4207875493, 4212110858, 
4215928942, 4219388367, 4222537267, 4225415529, 
4228056489, 4230488221, 4232734536, 4234815761, 
4236749356, 4238550398, 4240231977, 4241805513, 
4243281010, 4244667275, 4245972086, 4247202343, 
4248364183, 4249463089, 4250503969, 4251491232, 
4252428852, 4253320414, 4254169169, 4254978064, 
4255749782, 4256486769, 4257191258, 4257865296, 
4258510757, 4259129366, 4259722710, 4260292251, 
4260839340, 4261365226, 4261871068, 4262357938, 
4262826834, 4263278682, 4263714347, 4264134635, 
4264540296, 4264932033, 4265310504, 4265676324, 
4266030071, 4266372285, 4266703477, 4267024124, 
4267334677, 4267635559, 4267927170, 4268209888, 
4268484069, 4268750048, 4269008145, 4269258660, 
4269501878, 4269738070, 4269967490, 4270190383, 
4270406979, 4270617497, 4270822145, 4271021120, 
4271214612, 4271402799, 4271585851, 4271763932, 
4271937197, 4272105792, 4272269859, 4272429531, 
4272584939, 4272736203, 4272883440, 4273026764, 
4273166279, 4273302089, 4273434291, 4273562977, 
4273688238, 4273810159, 4273928821, 4274044303, 
4274156678, 4274266019, 4274372394, 4274475869, 
4274576505, 4274674363, 4274769501, 4274861972, 
4274951829, 4275039123, 4275123900, 4275206207, 
4275286087, 4275363581, 4275438730, 4275511570, 
4275582138, 4275650468, 4275716591, 4275780540, 
4275842342, 4275902026, 4275959617, 4276015139, 
4276068617, 4276120070, 4276169520, 4276216984, 
4276262481, 4276306026, 4276347634, 4276387317, 
4276425089, 4276460960, 4276494939, 4276527035, 
4276557254, 4276585601, 4276612082, 4276636699, 
4276659454, 4276680347, 4276699378, 4276716545, 
4276731843, 4276745269, 4276756815, 4276766475, 
4276774239, 4276780097, 4276784038, 4276786047, 
4276786109, 4276784209, 4276780328, 4276774446, 
4276766542, 4276756592, 4276744570, 4276730451, 
4276714204, 4276695798, 4276675201, 4276652375, 
4276627284, 4276599887, 4276570140, 4276537998, 
4276503413, 4276466332, 4276426702, 4276384463, 
4276339555, 4276291913, 4276241467, 4276188144, 
4276131868, 4276072556, 4276010120, 4275944470, 
4275875508, 4275803129, 4275727225, 4275647679, 
4275564367, 4275477159, 4275385915, 4275290486, 
4275190717, 4275086438, 4274977471, 4274863626, 
4274744700, 4274620477, 4274490724, 4274355195, 
4274213623, 4274065726, 4273911198, 4273749712, 
4273580915, 4273404430, 4273219848, 4273026727, 
4272824590, 4272612923, 4272391163, 4272158705, 
4271914885, 4271658984, 4271390211, 4271107707, 
4270810524, 4270497623, 4270167859, 4269819968, 
4269452547, 4269064041, 4268652718, 4268216642, 
4267753642, 4267261279, 4266736796, 4266177069, 
4265578544, 4264937154, 4264248226, 4263506362, 
4262705290, 4261837680, 4260894904, 4259866734, 
4258740950, 4257502823, 4256134423, 4254613696, 
4252913182, 4250998234, 4248824459, 4246333994, 
4243449892, 4240067428, 4236040118, 4231156253, 
4225097395, 4217360117, 4207095967, 4192747394, 
4171088088, 4134068812, 4053578179, 0
};

static const double zfn[ZIGLEVELS] = {
1, 0.97710143059782, 0.95987861830851, 0.94519830763172, 
0.93205927637968, 0.91999056461098, 0.90872536850199, 0.89809472692918, 
0.88798334868776, 0.87830823196574, 0.86900715699489, 0.86003198695311, 
0.85134452455213, 0.84291382270764, 0.8347143689422, 0.82672481886199, 
0.81892708786231, 0.81130568410951, 0.80384720853967, 0.79653997325589, 
0.7893737056314, 0.78233931561091, 0.77542871038842, 0.76863464512957, 
0.76195060148557, 0.75537068779605, 0.74888955640602, 0.74250233462364, 
0.73620456665171, 0.7299921644229, 0.72386136571625, 0.71780869827208, 
0.71183094888197, 0.70592513663134, 0.70008848962889, 0.69431842467991, 
0.6886125294583, 0.68296854680955, 0.67738436087973, 0.67185798481587, 
0.6663875498242, 0.66097129540631, 0.65560756062102, 0.65029477624233, 
0.64503145770312, 0.63981619872978, 0.63464766558632, 0.6295245918578, 
0.62444577371207, 0.619410065587, 0.61441637625712, 0.60946366523946, 
0.60455093950318, 0.59967725045204, 0.59484169115244, 0.59004339378279, 
0.58528152728292, 0.58055529518449, 0.57586393360571, 0.57120670939507, 
0.56658291841091, 0.56199188392465, 0.55743295513706, 0.55290550579771, 
0.54840893291904, 0.54394265557723, 0.53950611379255, 0.53509876748314, 
0.53072009548603, 0.52636959464053, 0.52204677892881, 0.51775117866965, 
0.51348233976116, 0.50923982296903, 0.50502320325688, 0.50083206915568, 
0.49666602216963, 0.49252467621576, 0.48840765709505, 0.48431460199287, 
0.48024515900682, 0.47619898670008, 0.47217575367868, 0.46817513819109, 
0.46419682774868, 0.46024051876584, 0.45630591621844, 0.45239273331948, 
0.44850069121096, 0.444629518671, 0.4407789518352, 0.43694873393151, 
0.43313861502772, 0.42934835179101, 0.42557770725868, 0.42182645061953, 
0.41809435700533, 0.41438120729177, 0.41068678790836, 0.40701089065686, 
0.40335331253773, 0.39971385558425, 0.39609232670375, 0.39248853752587, 
0.38890230425711, 0.38533344754175, 0.38178179232853, 0.37824716774291, 
0.37472940696465, 0.37122834711051, 0.3677438291216, 0.3642756976555, 
0.36082380098257, 0.35738799088659, 0.35396812256921, 0.35056405455838, 
0.34717564862021, 0.3438027696745, 0.34044528571341, 0.33710306772341, 
0.33377598961025, 0.33046392812682, 0.3271667628038, 0.32388437588306, 
0.32061665225349, 0.31736347938947, 0.31412474729155, 0.31090034842955, 
0.30769017768774, 0.30449413231222, 0.30131211186028, 0.29814401815175, 
0.29498975522227, 0.29184922927829, 0.28872234865397, 0.2856090237697, 
0.2825091670923, 0.27942269309683, 0.27634951822993, 0.27328956087471, 
0.27024274131705, 0.26720898171339, 0.26418820605984, 0.2611803401627, 
0.25818531161026, 0.25520304974587, 0.25223348564236, 0.24927655207758, 
0.2463321835112, 0.24340031606275, 0.24048088749069, 0.23757383717278, 
0.2346791060875, 0.23179663679657, 0.22892637342863, 0.22606826166403, 
0.22322224872061, 0.22038828334068, 0.21756631577903, 0.21475629779197, 
0.21195818262754, 0.20917192501673, 0.20639748116582, 0.20363480874977, 
0.2008838669068, 0.19814461623394, 0.19541701878386, 0.19270103806273, 
0.18999663902929, 0.18730378809511, 0.18462245312601, 0.18195260344482, 
0.17929420983526, 0.17664724454729, 0.17401168130367, 0.171387495308, 
0.16877466325412, 0.16617316333711, 0.16358297526567, 0.16100408027628, 
0.15843646114894, 0.1558801022247, 0.15333498942496, 0.15080111027283, 
0.14827845391637, 0.14576701115404, 0.14326677446239, 0.14077773802608, 
0.13829989777048, 0.13583325139685, 0.13337779842039, 0.13093354021134, 
0.12850048003918, 0.12607862312034, 0.12366797666949, 0.12126854995483, 
0.11888035435751, 0.11650340343556, 0.11413771299272, 0.11178330115245, 
0.10944018843764, 0.10710839785639, 0.10478795499442, 0.1024788881147, 
0.10018122826488, 0.097895009393179, 0.095620268473705, 0.093357045641809, 
0.091105384340662, 0.088865331480041, 0.086636937608583, 0.084420257100895, 
0.082215348361065, 0.080022274044367, 0.077841101299147, 0.075671902031198, 
0.073514753193218, 0.071369737102363, 0.069236941789315, 0.06711646138285, 
0.065008396534479, 0.062912854888512, 0.060829951603757, 0.058759809934132, 
0.056702561876794, 0.054658348897885, 0.052627322747973, 0.050609646381549, 
0.048605494997894, 0.046615057224181, 0.044638536466257, 0.042676152458222, 
0.040728143049248, 0.038794766275447, 0.036876302776783, 0.034973058635058, 
0.033085368730228, 0.031213600740953, 0.029358159954282, 0.027519495103399, 
0.02569810552843, 0.023894550064357, 0.022109458219765, 0.020343544449317, 
0.018597626690655, 0.016872650919093, 0.015169724428837, 0.013490162179794, 
0.011835553465957, 0.010207861685007, 0.008609581203028, 0.0070440001618338, 
0.0055156798328771, 0.0040314406569058, 0.0026028041937557, 0.0012544610762767
};

static const double zwn[ZIGLEVELS] = {
9.1478851695481e-10, 5.0115208832011e-11, 6.6630615736164e-11, 7.8170456525769e-11, 
8.7340273839945e-11, 9.5086087855792e-11, 1.0186831870851e-10, 1.079489492412e-10, 
1.1349259815484e-10, 1.186101026495e-10, 1.2337999544575e-10, 1.2786014996961e-10, 
1.3209456738668e-10, 1.3611756277377e-10, 1.3995646787861e-10, 1.4363344318803e-10, 
1.4716673191116e-10, 1.5057155148405e-10, 1.5386074229221e-10, 1.5704524939719e-10, 
1.6013448669183e-10, 1.6313661656548e-10, 1.6605876773433e-10, 1.6890720707322e-10, 
1.716874767228e-10, 1.7440450463239e-10, 1.7706269453364e-10, 1.7966599981084e-10, 
1.8221798463526e-10, 1.8472187493298e-10, 1.8718060116695e-10, 1.8959683447525e-10, 
1.919730173773e-10, 1.943113900078e-10, 1.9661401264493e-10, 1.9888278514931e-10, 
2.0111946381326e-10, 2.0332567602737e-10, 2.0550293309828e-10, 2.0765264149312e-10, 
2.097761127391e-10, 2.118745721685e-10, 2.1394916666864e-10, 2.1600097157111e-10, 
2.1803099679343e-10, 2.2004019232975e-10, 2.2202945317215e-10, 2.2399962373315e-10, 
2.2595150182919e-10, 2.2788584227718e-10, 2.2980336014879e-10, 2.3170473372137e-10, 
2.3359060715914e-10, 2.3546159295419e-10, 2.3731827415288e-10, 2.3916120639031e-10, 
2.4099091975259e-10, 2.4280792048441e-10, 2.4461269255729e-10, 2.4640569911216e-10, 
2.4818738378843e-10, 2.4995817195017e-10, 2.5171847181915e-10, 2.5346867552316e-10, 
2.5520916006735e-10, 2.5694028823533e-10, 2.5866240942634e-10, 2.6037586043383e-10, 
2.620809661706e-10, 2.6377804034486e-10, 2.6546738609138e-10, 2.6714929656126e-10, 
2.6882405547385e-10, 2.7049193763359e-10, 2.7215320941476e-10, 2.738081292165e-10, 
2.7545694789034e-10, 2.7709990914247e-10, 2.7873724991253e-10, 2.803692007307e-10, 
2.8199598605465e-10, 2.8361782458787e-10, 2.852349295807e-10, 2.8684750911524e-10, 
2.8845576637534e-10, 2.9005989990275e-10, 2.9166010384025e-10, 2.9325656816283e-10, 
2.9484947889764e-10, 2.9643901833348e-10, 2.9802536522063e-10, 2.9960869496149e-10, 
3.011891797929e-10, 3.027669889605e-10, 3.0434228888573e-10, 3.0591524332598e-10, 
3.074860135284e-10, 3.0905475837765e-10, 3.1062163453817e-10, 3.1218679659121e-10, 
3.1375039716712e-10, 3.1531258707309e-10, 3.1687351541674e-10, 3.1843332972583e-10, 
3.1999217606444e-10, 3.2155019914576e-10, 3.2310754244186e-10, 3.246643482906e-10, 
3.2622075800001e-10, 3.2777691195019e-10, 3.2933294969319e-10, 3.3088901005074e-10, 
3.3244523121042e-10, 3.3400175082003e-10, 3.3555870608073e-10, 3.3711623383883e-10, 
3.3867447067658e-10, 3.4023355300202e-10, 3.4179361713812e-10, 3.4335479941134e-10, 
3.449172362397e-10, 3.4648106422067e-10, 3.4804642021884e-10, 3.4961344145369e-10, 
3.511822655875e-10, 3.5275303081355e-10, 3.5432587594486e-10, 3.5590094050346e-10, 
3.5747836481053e-10, 3.5905829007737e-10, 3.6064085849751e-10, 3.6222621334e-10, 
3.6381449904416e-10, 3.6540586131581e-10, 3.670004472253e-10, 3.6859840530739e-10, 
3.7019988566324e-10, 3.718050400646e-10, 3.7341402206059e-10, 3.7502698708696e-10, 
3.7664409257835e-10, 3.7826549808356e-10, 3.7989136538412e-10, 3.8152185861644e-10, 
3.831571443977e-10, 3.847973919558e-10, 3.8644277326365e-10, 3.8809346317801e-10, 
3.897496395833e-10, 3.914114835406e-10, 3.9307917944223e-10, 3.9475291517225e-10, 
3.9643288227326e-10, 3.9811927611991e-10, 3.9981229609958e-10, 4.0151214580058e-10, 
4.0321903320851e-10, 4.0493317091112e-10, 4.0665477631241e-10, 4.083840718564e-10, 
4.1012128526127e-10, 4.1186664976452e-10, 4.1362040437995e-10, 4.1538279416707e-10, 
4.1715407051395e-10, 4.1893449143423e-10, 4.2072432187938e-10, 4.225238340672e-10, 
4.2433330782756e-10, 4.261530309668e-10, 4.2798329965182e-10, 4.2982441881539e-10, 
4.3167670258426e-10, 4.3354047473143e-10, 4.3541606915463e-10, 4.3730383038276e-10, 
4.3920411411232e-10, 4.4111728777626e-10, 4.4304373114743e-10, 4.449838369796e-10, 
4.4693801168865e-10, 4.4890667607729e-10, 4.5089026610667e-10, 4.5288923371877e-10, 
4.5490404771345e-10, 4.5693519468507e-10, 4.5898318002327e-10, 4.6104852898369e-10, 
4.6313178783442e-10, 4.6523352508497e-10, 4.6735433280495e-10, 4.6949482804066e-10, 
4.7165565433842e-10, 4.7383748338474e-10, 4.7604101677407e-10, 4.7826698791667e-10, 
4.805161641e-10, 4.8278934871908e-10, 4.8508738369253e-10, 4.8741115208372e-10, 
4.8976158094802e-10, 4.9213964443049e-10, 4.9454636714088e-10, 4.9698282783659e-10, 
4.9945016344811e-10, 5.0194957348615e-10, 5.0448232487504e-10, 5.0704975726297e-10, 
5.096532888673e-10, 5.1229442292105e-10, 5.1497475479687e-10, 5.176959798963e-10, 
5.2045990240549e-10, 5.2326844503496e-10, 5.2612365987966e-10, 5.2902774055865e-10, 
5.3198303582032e-10, 5.349920648321e-10, 5.3805753441236e-10, 5.4118235850975e-10, 
5.4436968029303e-10, 5.4762289728469e-10, 5.5094569005839e-10, 5.5434205512748e-10, 
5.5781634278442e-10, 5.6137330081797e-10, 5.6501812524383e-10, 5.6875651945056e-10, 
5.7259476350115e-10, 5.7653979576767e-10, 5.8059930964246e-10, 5.8478186881163e-10, 
5.8909704555766e-10, 5.9355558786838e-10, 5.981696229009e-10, 6.0295290676928e-10, 
6.0792113397877e-10, 6.1309232454029e-10, 6.1848731352045e-10, 6.2413037753653e-10, 
6.3005004712642e-10, 6.3628017568877e-10, 6.4286136930223e-10, 6.4984293499996e-10, 
6.5728559197658e-10, 6.652653367448e-10, 6.7387910993869e-10, 6.8325338236121e-10, 
6.9355768403019e-10, 7.0502695963406e-10, 7.1800075173536e-10, 7.3299724341921e-10, 
7.5086784216768e-10, 7.7316823738581e-10, 8.0326004350637e-10, 8.5109388898754e-10
};

/* BEGIN MT19937 CODE */

/* seed the mt19937 random number generator */
void mt19937_srand(it_uint32_t seed)
{
  state[0]= seed & 0xffffffffUL;
  for(idx = 1; idx < N; idx++)
    {
      state[idx] = (1812433253UL * (state[idx-1] ^ (state[idx-1] >> 30)) + idx); 
      /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
      /* In the previous versions, MSBs of the seed affect   */
      /* only MSBs of the array state[].                        */
      /* 2002/01/09 modified by Makoto Matsumoto             */
      state[idx] &= 0xffffffffUL;  /* for >32 bit machines */
    }

  left = 1; 
  initf = 1; 

  return;
}

void mt19937_srand_by_array(it_uint32_t init_key[], it_uint32_t key_length)
{
    int i, j, k;

    mt19937_srand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        state[i] = (state[i] ^ ((state[i-1] ^ (state[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        state[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { state[0] = state[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        state[i] = (state[i] ^ ((state[i-1] ^ (state[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        state[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { state[0] = state[N-1]; i=1; }
    }

    state[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
    left = 1; initf = 1;

    return;
}

static void mt19937_next_state(void)
{
  it_uint32_t *p=state;
  int j;

  /* if init_genrand() has not been called, */
  /* a default initial seed is used         */
  if (initf==0) mt19937_srand(5489UL);
  
  left = N;
  next = state;
  
  for (j=N-M+1; --j; p++) 
    *p = p[M] ^ TWIST(p[0], p[1]);
  
  for (j=M; --j; p++) 
    *p = p[M-N] ^ TWIST(p[0], p[1]);
  
  *p = p[M-N] ^ TWIST(p[0], state[0]);
  
  return;
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned int mt19937_rand_int32(void)
{
  it_uint32_t y;
  
  if (--left == 0) mt19937_next_state();
  y = *next++;
  
  /* Tempering */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);
  
  return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
int mt19937_rand_int31(void)
{
  it_uint32_t y;

  if (--left == 0) mt19937_next_state();
  y = *next++;
  
  /* Tempering */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);
  
  return (long)(y>>1);
}

/* generates a random number on [0,1]-real-interval */
double mt19937_rand_real1(void)
{
  it_uint32_t y;

  if (--left == 0) mt19937_next_state();
  y = *next++;
  
  /* Tempering */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);
  
  return (double)y * (1.0/4294967295.0); 
  /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double mt19937_rand_real2(void)
{
  it_uint32_t y;

  if (--left == 0) mt19937_next_state();
  y = *next++;
  
  /* Tempering */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);
  
  return (double)y * (1.0/4294967296.0); 
  /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double mt19937_rand_real3(void)
{
  it_uint32_t y;

  if (--left == 0) mt19937_next_state();
  y = *next++;
  
  /* Tempering */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);
  
  return ((double)y + 0.5) * (1.0/4294967296.0); 
  /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double mt19937_rand_res53(void) 
{ 
  it_uint32_t a=mt19937_rand_int32()>>5, b=mt19937_rand_int32()>>6; 
  return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */

/* -------------------------------------------------------------------------- */

/* generate a random number using the mt19937 algorithm */
/* this number has a precision of 32 bits and is stored */
/* in double-precision in the [0, 1) range.             */
static double mt19937_rand()
{
  
  return( mt19937_rand_real2() ); 
  
}

/* END MT19937 CODE */

/* create and initialize the random number generator */
void it_randomize(void)
{
  unsigned int usec;

#ifdef WIN32
  /* the windows version */
  SYSTEMTIME SystemTime;

  GetLocalTime(&SystemTime);
  usec = (unsigned int) SystemTime.wMilliseconds;
#else
  /* the unix version */
  struct timeval tp;

  /* seed the random number generator with the current time */
  gettimeofday(&tp, NULL);
  usec = tp.tv_usec;
#endif

  it_seed(usec);
}

void it_seed(int seed)
{
  mt19937_srand((it_uint32_t) seed);
}

/* generate a random value uniformly distributed in [0,1) */
double it_rand(void)
{
  return(mt19937_rand());
}

/*
  Use this method to generate Gaussian distributed values 
  with N(0,1) model. This code closely follows original 
  implementation from Marsaglia and Tsang. 
 */
double it_randn( void ) 
{

  unsigned long int i, j; 
  double x, y; 
  int s; 

  while ( 1 ) 
    {

      j = mt19937_rand_int32(); 

      i = j & 0x000000FF; /* "Form i from the last 8 bits of j" */
      
      x = j*zwn[i]; 

      s = j & 0x00000800 ? 1. : -1.; /* Get a sign from any bit of j */

      if ( j < zkn[i] ) /* Should exit here 99% of the time */ 
	return( s*x ); 

      
      if ( !i ) /* "return an x from the tail" */ 
	{
	  do {
	    x = -log( mt19937_rand_real3() ) * ZIGRINV; 
	    y = -log( mt19937_rand_real3() ); 
	  } while( y+y < x*x ); 

	  return( s*(ZIGR+x) );
	}

      if ( mt19937_rand_real3()*(zfn[i-1]-zfn[i]) < exp(-.5*x*x)-zfn[i] ) 
	return( s*x ); 
    }

}

/* generate a random variable from its probability
   density function using the acceptance-rejection method.
   the pdf is assumed to be zero outside [a, b]
   and centered on its maximum value.
*/
double it_randpdf(double a, double b, it_function_t pdf, it_args_t args)
{
  double x, y;
  double max;

  max = pdf(0, args); /* evaluate the maximum value of the pdf. */
  do {
    /* draw x and y uniformly in [a,b]x[0,1] */
    x = a + it_rand() * (b - a);
    y = max * it_rand();

    /* check if (x,y) is below the pdf curve */
    if(y < pdf(x, args)) return(x);
  } while(1);
}



int it_rand_memoryless( vec pdf ) 
{
  double x;
  int s, e, h, r; /* r is the result */
  vec pdfcs = vec_cum_sum( pdf );
  vec_ins( pdfcs, 0, 0.0 );

  /* generate a uniformly distributed value  */
  /* and search the integral of the pdf for  */
  /* the abciss corresponding to that random */
  /* value.                                  */
  x = it_rand();
  s = 0;
  e = vec_length(pdfcs) - 1;

  while(e > s + 1) {
    /* compare with the median of the range */
    /* and narrow down the range. */
    h = (s + e) / 2; 

    if(x < pdfcs[h])
      e = h;
    else
      s = h;
  }

  if(x < pdfcs[e])
    r = s;
  else
    r = e;

  vec_delete( pdfcs );
  return r;
}
