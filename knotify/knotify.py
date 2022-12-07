#
# Copyright © 2022 Christos Pavlatos, George Rassias, Christos Andrikos,
#                  Evangelos Makris, Aggelos Kolaitis
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the “Software”), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

from typing import Tuple, List

from oslo_config import cfg, types

from knotify import hairpin
from knotify.algorithm.base import BaseAlgorithm
from knotify.algorithm.ipknot import IPKnot
from knotify.algorithm.knotify import Knotify
from knotify.algorithm.knotty import Knotty
from knotify.algorithm.ihfold import IHFold
from knotify.algorithm.hotknots import HotKnots
from knotify.energy.vienna import ViennaEnergy
from knotify.energy.external import ExternalEnergy
from knotify.energy.pkenergy import PKEnergy
from knotify.pairalign.cpairalign import CPairAlign
from knotify.pairalign.bulges import BulgesPairAlign
from knotify.extensions.skip_final_au import SkipFinalAU
from knotify.parsers.yaep import YaepParser
from knotify.parsers.bruteforce import BruteForceParser

CSV_OPTS = [
    cfg.StrOpt("csv"),
    cfg.StrOpt("results-csv"),
]

PSEUDOKNOT_OPTS = [
    cfg.StrOpt("parser", choices=["bruteforce", "yaep"], default="yaep"),
    cfg.StrOpt("yaep-library-path", default="./libpseudoknot.so"),
    cfg.StrOpt("bruteforce-library-path", default="./libbruteforce.so"),
    cfg.BoolOpt("allow-ug", default=False),
    cfg.IntOpt("max-dd-size", default=2),
    cfg.IntOpt("min-dd-size", default=0),
    cfg.IntOpt("max-window-size", default=204),
    cfg.IntOpt("min-window-size", default=6),
    cfg.FloatOpt("max-window-size-ratio", default=0),
    cfg.FloatOpt("min-window-size-ratio", default=0),
    cfg.IntOpt("max-stem-allow-smaller", default=2),
    cfg.BoolOpt("prune-early", default=True),
]

PAIRALIGN_OPTS = [
    cfg.ListOpt(
        "pairalign",
        item_type=types.String(choices=["consecutive", "bulges"]),
        default=["consecutive"],
    ),
    cfg.BoolOpt("allow-skip-final-au", default=False),
    cfg.StrOpt("skip-final-au-library-path", default="./libskipfinalau.so"),
    cfg.StrOpt("consecutive-pairalign-library-path", default="./libcpairalign.so"),
    cfg.StrOpt("bulges-library-path", default="./libbulges.so"),
    cfg.IntOpt("max-bulge-size", default=1),
    cfg.IntOpt("min-stems-after-bulge", default=1),
    cfg.BoolOpt("symmetric-bulges", default=True),
    cfg.BoolOpt("count-stems-from-bulges", default=False),
]

HAIRPIN_OPTS = [
    cfg.StrOpt("hairpin-grammar"),
    cfg.BoolOpt("hairpin-allow-ug"),
    cfg.IntOpt("min-hairpin-size", default=hairpin.MIN_HAIRPIN_SIZE),
    cfg.IntOpt("min-hairpin-stems", default=hairpin.MIN_HAIRPIN_STEMS),
    cfg.IntOpt("max-hairpins-per-loop", default=hairpin.MAX_HAIRPINS_PER_LOOP),
    cfg.IntOpt("max-hairpin-bulge", default=hairpin.MAX_HAIRPIN_BULGE),
]

ENERGY_OPTS = [
    cfg.StrOpt("energy", default="vienna", choices=["vienna", "pkenergy", "external"]),
    cfg.StrOpt("pkenergy", default="./libpkenergy.so"),
    cfg.StrOpt("pkenergy-config-dir", default="./pkenergy/hotknots/params"),
    cfg.StrOpt(
        "pkenergy-model",
        default="dp",
        choices=["dp", "re", "cc2006a", "cc2006b", "cc2006c"],
    ),
    cfg.StrOpt("external-energy-executable"),
]

ALGORITHM_OPTS = [
    cfg.StrOpt(
        "algorithm",
        default="knotify",
        choices=["knotify", "ipknot", "knotty", "hotknots", "ihfold"],
    ),
    cfg.StrOpt("ipknot-executable", default="./.ipknot/ipknot"),
    cfg.StrOpt("knotty-executable", default="./.knotty/knotty"),
    cfg.StrOpt("ihfold-executable", default="./.iterative-hfold/HFold_iterative"),
    cfg.StrOpt("hotknots-dir", default="./.hotknots/HotKnots_v2.0"),
]


class ConfigOpts(cfg.ConfigOpts):
    """
    A cfg.ConfigOpts with added type-hints for the available configuration options.

    Configuration can also be loaded from environment variables, using KNOTIFY_$OPTION_NAME,
    e.g. KNOTIFY_YAEP_LIBRARY_PATH.
    """

    def __init__(self, *args, **kwargs):
        super(ConfigOpts, self).__init__(*args, **kwargs)
        self._env_driver.get_name = lambda _, opt: "KNOTIFY_{}".format(opt.upper())

    # PSEUDOKNOT_OPTS
    parser: str
    yaep_library_path: str
    bruteforce_library_path: str
    allow_ug: bool
    max_dd_size: int
    min_dd_size: int
    max_window_size: int
    min_window_size: int
    max_window_size_ratio: float
    min_window_size_ratio: float
    max_stem_allow_smaller: int
    prune_early: bool

    # PAIRALIGN_OPTS
    pairalign: List[str]
    allow_skip_final_au: bool
    skip_final_au_library_path: str
    consecutive_pairalign_library_path: str
    bulges_library_path: str
    max_bulge_size: int
    min_stems_after_bulge: int
    symmetric_bulges: bool
    count_stems_from_bulges: bool

    # HAIRPIN_OPTS
    hairpin_grammar: str
    hairpin_allow_ug: bool
    min_hairpin_size: int
    min_hairpin_stems: int
    max_hairpins_per_loop: int
    max_hairpin_bulge: int

    # ENERGY_OPTS
    energy: str
    pkenergy: str
    pkenergy_config_dir: str
    pkenergy_model: str
    external_energy_executable: str

    # ALGORITHM_OPTS
    algorithm: str
    ipknot_executable: str
    knotty_executable: str
    ihfold_executable: str
    hotknots_dir: str


def new_options() -> ConfigOpts:
    """
    Construct a new options object.
    """
    c = ConfigOpts()

    for options in [
        CSV_OPTS,
        PSEUDOKNOT_OPTS,
        PAIRALIGN_OPTS,
        ALGORITHM_OPTS,
        ENERGY_OPTS,
        HAIRPIN_OPTS,
    ]:
        c.register_opts(options)
        c.register_cli_opts(options)

    return c


def from_options(opts: ConfigOpts) -> Tuple[BaseAlgorithm, dict]:
    """
    Initialize the knotify execution engine and a dict with configuration.
    """
    if opts.hairpin_allow_ug is None:
        opts.hairpin_allow_ug = opts.allow_ug

    rna_parser_args = {
        "max_dd_size": opts.max_dd_size,
        "min_dd_size": opts.min_dd_size,
        "max_window_size": opts.max_window_size,
        "min_window_size": opts.min_window_size,
        "max_window_size_ratio": opts.max_window_size_ratio,
        "min_window_size_ratio": opts.min_window_size_ratio,
        "allow_ug": opts.allow_ug,
    }
    parser = None
    if opts.parser == "yaep":
        parser = YaepParser(opts.yaep_library_path, **rna_parser_args)
    elif opts.parser == "bruteforce":
        parser = BruteForceParser(opts.bruteforce_library_path, **rna_parser_args)

    energy = None
    if opts.energy == "vienna":
        energy = ViennaEnergy()
    elif opts.energy == "pkenergy":
        energy = PKEnergy(opts.pkenergy, opts.pkenergy_config_dir, opts.pkenergy_model)
    elif opts.energy == "external":
        energy = ExternalEnergy(opts.external_energy_executable)

    algorithm = None
    if opts.algorithm == "knotify":
        algorithm = Knotify()
    elif opts.algorithm == "ipknot":
        algorithm = IPKnot()
    elif opts.algorithm == "knotty":
        algorithm = Knotty()
    elif opts.algorithm == "ihfold":
        algorithm = IHFold()
    elif opts.algorithm == "hotknots":
        algorithm = HotKnots()

    pairalign = []
    if "consecutive" in opts.pairalign:
        pairalign.append(CPairAlign(opts.consecutive_pairalign_library_path).pairalign)
    if "bulges" in opts.pairalign:
        pairalign.append(
            BulgesPairAlign(
                opts.max_bulge_size,
                opts.min_stems_after_bulge,
                opts.symmetric_bulges,
                opts.count_stems_from_bulges,
                library_path=opts.bulges_library_path,
            ).pairalign
        )

    skip_final_au = None
    if opts.allow_skip_final_au:
        skip_final_au = SkipFinalAU(opts.skip_final_au_library_path).get_candidates

    return algorithm, {
        "parser": parser.detect_pseudoknots,
        "csv": opts.csv,
        "allow_ug": opts.allow_ug,
        "pairalign": pairalign,
        "skip_final_au": skip_final_au,
        "allow_skip_final_au": opts.allow_skip_final_au,
        "max_stem_allow_smaller": opts.max_stem_allow_smaller,
        "prune_early": opts.prune_early,
        "hairpin_grammar": opts.hairpin_grammar,
        "hairpin_allow_ug": opts.hairpin_allow_ug,
        "min_hairpin_size": opts.min_hairpin_size,
        "min_hairpin_stems": opts.min_hairpin_stems,
        "max_hairpin_bulge": opts.max_hairpin_bulge,
        "max_hairpins_per_loop": opts.max_hairpins_per_loop,
        "energy": energy,
        "ipknot_executable": opts.ipknot_executable,
        "knotty_executable": opts.knotty_executable,
        "ihfold_executable": opts.ihfold_executable,
        "hotknots_dir": opts.hotknots_dir,
    }
