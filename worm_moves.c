// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#include "worm_moves.h"

// ----------------------------------------------------------------
char* worm_move_to_desc(int worm_move)
{
	switch (worm_move) {
	case     WORM_OPEN:      return "open";      break;
	case     WORM_CLOSE:     return "close";     break;
	case     WORM_HEAD_SWAP: return "head swap"; break;
	case     WORM_TAIL_SWAP: return "tail swap"; break;
	default:                 return "???";       break;
	}
}
