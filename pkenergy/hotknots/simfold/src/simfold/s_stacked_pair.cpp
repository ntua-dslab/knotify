/***************************************************************************
                          s_stacked_pair.cpp  -  description
                             -------------------
    begin                : Thu Apr 11 2002
    copyright            : (C) 2002 by Mirela Andronescu
    email                : andrones@cs.ubc.ca
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <string.h>

#include "common.h"
#include "constants.h"
#include "externs.h"
#include "params.h"
#include "s_stacked_pair.h"
#include "simfold.h"
#include "structs.h"

s_stacked_pair::s_stacked_pair(int *seq, int length)
// The constructor
{
  sequence = seq; // just refer it from where it is in memory
  seqlen = length;
  this->V = NULL;
}

s_stacked_pair::~s_stacked_pair()
// The destructor
{}

PARAMTYPE s_stacked_pair::compute_energy(int i, int j)
// compute the free energy of the structure closed by this stacked pair
{
  PARAMTYPE min = INF, local_energy, V_energy;

  V_energy = V->get_energy(i + 1, j - 1);

  local_energy =
      stack[sequence[i]][sequence[j]][sequence[i + 1]][sequence[j - 1]];

  // TODO
  //     if (sequence[i] == 0 && sequence[j] == 3 && sequence[i+1] == 0 &&
  //     sequence[j-1] == 3)
  //         local_energy = stack[0][3][0][3];
  //     else if (sequence[j-1] == 0 && sequence[i+1] == 3 && sequence[j] == 0
  //     && sequence[i] == 3)
  //         local_energy = stack[0][3][0][3];
  //     else
  //         local_energy = -100;

  min = V_energy + local_energy;
  return min;
}

PARAMTYPE s_stacked_pair::get_energy(int i, int j, int *sequence)
// returns the free energy of the stacked pair closed at (i,j)
{
  if (i + 1 >= j - 1)
    return INF;
  // TODO
  //     if (sequence[i] == 0 && sequence[j] == 3 && sequence[i+1] == 0 &&
  //     sequence[j-1] == 3)
  //         return stack[0][3][0][3];
  //     else if (sequence[j-1] == 0 && sequence[i+1] == 3 && sequence[j] == 0
  //     && sequence[i] == 3)
  //         return stack[0][3][0][3];
  //     else
  //         return -100;

  return stack[sequence[i]][sequence[j]][sequence[i + 1]][sequence[j - 1]];
}

PARAMTYPE s_stacked_pair::get_enthalpy(int i, int j, int *sequence)
// returns the enthalpy of the stacked pair closed at (i,j)
{
  if (i + 1 >= j - 1)
    return INF;
  return enthalpy_stack[sequence[i]][sequence[j]][sequence[i + 1]]
                       [sequence[j - 1]];
}

void s_stacked_pair::count_get_energy(int i, int j, int *sequence,
                                      double *counter)
// used for parameter learning, not for folding
// update the counter vector accordingly
// Mirela: Nov 23, 2003
{
  char type[100];
  int index;
  if (1000 * sequence[i] + 100 * sequence[j] + 10 * sequence[i + 1] +
          sequence[j - 1] >
      1000 * sequence[j - 1] + 100 * sequence[i + 1] + 10 * sequence[j] +
          sequence[i])
    sprintf(type, "stack[%d][%d][%d][%d]", sequence[j - 1], sequence[i + 1],
            sequence[j], sequence[i]);
  else
    sprintf(type, "stack[%d][%d][%d][%d]", sequence[i], sequence[j],
            sequence[i + 1], sequence[j - 1]);
  index = structure_type_index(type);
  counter[index]++;
  return;
}
