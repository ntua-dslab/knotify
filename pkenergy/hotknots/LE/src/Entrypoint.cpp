// Author: Angelos Kolaitis
// Email: neoaggelos@gmail.com
// Date: 2021-10-14
// Description: Entrypoint for the PK energy calculation algorithm. Defined
//              with C linkage, to allow loading from foreign code.

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <stdio.h>
#include <string.h>

#include "Bands.h"
#include "Input.h"
#include "Loop.h"
#include "Stack.h"
#include "commonPK.h" // detect_original_PKed_pairs_many()
#include "init.h"     // init_data()
#include "initPK.h"   // init_dataPK()
#include "params.h"   // fill_data_structures_with_new_parameters(), RNA

extern "C" {

// initialize is called once before starting energy calculations to load
// parameters. It should point to the `hotknots/params` directory of this
// repository.
void initialize(char *config_dir) {
  char multirnafold[200], pkenergy[200], constrdangles[200];
  snprintf(multirnafold, 200, "%s/multirnafold.conf", config_dir);
  snprintf(pkenergy, 200, "%s/pkenergy.conf", config_dir);
  snprintf(constrdangles, 200, "%s/turner_parameters_fm363_constrdangles.txt",
           config_dir);

  init_data("", multirnafold, RNA, 37);
  fill_data_structures_with_new_parameters(constrdangles);
  init_dataPK("", pkenergy, RNA, 37);
}

// Given an RNA sequence and its secondary structure, calculate the MFE
double get_energy(char *sequence, char *structure) {
  int size = strlen(structure);
  short pairseq[size + 1];
  detect_original_PKed_pairs_many(structure, pairseq);

  ReadInput *R = new ReadInput(size, sequence, pairseq);
  Stack *s = new Stack(R);
  Bands *B = new Bands(R, s);
  Loop *L = new Loop(0, MaxN + 1, R, B, s);

  for (int i = 1; i <= R->Size; i++) {
    if (R->BasePair(i) >= 0) {
      int a, b; // will storethe borders of a closed regoin
      if (s->Add(i, a, b)) {
        // If a closed region is identifed add it to the tree by calling addLoop
        L->addLoop(a, b);
      };
    };
  }
  double result =
      (-L->EnergyViaSimfold(DP) - L->EnergyDanglingViaSimfold(DP)) / 1000;

  delete R, s, B, L;
  return result;
}
}
