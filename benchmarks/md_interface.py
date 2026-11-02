"""Replicated interface-depth test (run F): seeds of many sizes and planar interfaces.

Purpose: measure the probability that a labeler assigns the correct bulk
polytype as a function of depth below a constructed ice/liquid interface, on
independent preparations, and test whether that curve predicts the labeled
nucleus size of seed sizes not used to fit it.

Systems, TIP4P/Ice at 220 K, three independent preparations p1-p3:
  F_seed_<kind>_<N>_p<k>   spherical Ih, Ic or stacking-disordered (Isd) seed of
                           N = 50, 100, 200, 400, 800 or 1600 molecules, cut
                           around a random lattice site of a large crystal with
                           a new proton configuration, randomly rotated, and
                           inserted at the centre of a 32,768-molecule liquid
  F_slab_<kind>_p<k>       planar control: a periodic Ih, Ic or Isd slab (basal
                           plane normal to z) in contact with liquid on both faces
Each preparation also has its own liquid (tiled seed-test liquid, remelted at
300 K with a new velocity seed and re-equilibrated at 220 K, NVT).

Confirmation set (run G, preparation p4, built and analyzed only after the
calibration was frozen in confirmation_protocol.md):
  G_seed_<kind>_<N>_p4     seeds of N = 300 and 600 molecules, sizes not used before
  G_prism_<kind>_p4        planar slab with primary prism faces {10-10} (stacking
                           axis in the interface plane), standard protocol
  G_long_<kind>_p4         basal slab, frames at 2-20 ps every 2 ps and then
                           every 20 ps to 500 ps (time dependence of the dead layer)

Protocol as the seed test (ice_nucleation_pg/data/seedtest/md.py): crystal
oxygens restrained (5000 kJ/mol/nm^2) during minimization and 20 ps, so the
liquid wets the crystal; released; 10 frames 2 ps apart. Construction labels
(seed membership, depth below the constructed surface, layer stacking) are
saved per molecule in <system>.truth.npz.

Needs OpenMM with CUDA (e.g. the gencompas conda environment):
  python md_interface.py --list
  python md_interface.py prepare            # structures (CPU) and liquids (GPU)
  python md_interface.py all | <systems>
  python md_interface.py prepare-confirmation && python md_interface.py confirmation
Outputs go to cache/interface/ (not committed).
"""

import argparse
import json
from pathlib import Path
import sys
import time

import networkx as nx
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

import md_replicas
from md_replicas import A, B, SHIFT, proton_configuration, water_sites, write_gro

SEEDTEST = md_replicas.SEEDTEST
OUT = Path(__file__).resolve().parent/"cache/interface"
sys.path.insert(0, str(SEEDTEST))
KINDS = ("Ih", "Ic", "Isd")
SIZES = (50, 100, 200, 400, 800, 1600)
PREPS = (1, 2, 3)
CRYSTAL = (12, 7, 18)  # NX, NY, layers: 54 x 55 x 66 A, large enough for R ~ 23 A
SLAB = (12, 7, 12)  # slab of 12 layers (44 A) under 42 A of liquid
LIQUID_TILE = 2  # seed-test liquid (4096 molecules, 50.5 A) tiled 2x2x2
INSERT_MARGIN = 2.2  # A; liquid removed within R + margin, as in the seed test
OVERLAP = 2.5  # A; minimum O-O distance when liquid is cut into a new periodic box
RESTRAINT, RESTRAINED_PS, FRAMES, STRIDE_PS = 5000., 20, 10, 2
LIQUID_MELT_PS, LIQUID_EQUIL_PS = 50, 300
CONFIRM_PREP, CONFIRM_SIZES = 4, (300, 600)
PRISM = (12, 6, 12)  # NX, NY, layers; the slab normal is y (46.8 A), the in-plane period 54 x 44 A
LONG_TIMES = (*range(2, 21, 2), *range(40, 501, 20))  # ps after release


def systems():
    out = [f"F_seed_{k}_{n}_p{p}" for p in PREPS for k in KINDS for n in SIZES]
    return out + [f"F_slab_{k}_p{p}" for p in PREPS for k in KINDS]


def confirmation_systems():
    p = CONFIRM_PREP
    return ([f"G_seed_{k}_{n}_p{p}" for k in KINDS for n in CONFIRM_SIZES]
            + [f"G_prism_{k}_p{p}" for k in KINDS] + [f"G_long_{k}_p{p}" for k in KINDS])


def stacking(kind, layers, rng):
    if kind == "Ih":
        return [0, 1]*(layers//2)
    if kind == "Ic":
        return [0, 1, 2]*(layers//3)
    while True:  # random, consecutive layers differ (including the wrap), about half hexagonal
        s = [0]
        for _ in range(layers-1):
            s.append(int(rng.choice([x for x in range(3) if x != s[-1]])))
        hexagonal = [s[k-1] == s[(k+1) % layers] for k in range(layers)]
        if s[-1] != s[0] and .4*layers <= sum(hexagonal) <= .6*layers:
            return s


def crystal(kind, shape, rng):
    """Proton-disordered crystal: O,H,H,M sites, box, per-molecule stacking (1 hex, 2 cubic), layer, bonds."""
    nx_, ny, layers = shape
    sequence = stacking(kind, layers, rng)
    box = np.array([nx_*A, ny*A*np.sqrt(3), layers*4*B/3])
    pos, layer = [], []
    for k in range(layers):
        for ix in range(nx_):
            for iy in range(ny):
                for bx, by in ((0., 0.), (A/2, A*np.sqrt(3)/2)):
                    for dz in (0., B):
                        pos.append((ix*A + bx, iy*A*np.sqrt(3) + by + SHIFT[sequence[k]], k*4*B/3 + dz))
                        layer.append(k)
    o = np.mod(np.array(pos), box)
    hexagonal = np.array([sequence[k-1] == sequence[(k+1) % layers] for k in range(layers)])
    stack = np.where(hexagonal[np.array(layer)], 1, 2).astype(np.int8)
    donors, dipole, _ = proton_configuration(o, box, rng)
    pairs = cKDTree(o, boxsize=box).query_pairs(B*1.1, output_type="ndarray")
    return water_sites(o, box, donors), box, stack, np.array(layer), pairs, "".join("ABC"[x] for x in sequence), dipole


def read_gro(path):
    lines = Path(path).read_text().splitlines()
    xyz = np.array([[float(l[20:28]), float(l[28:36]), float(l[36:44])] for l in lines[2:-1]])*10
    return xyz.reshape(-1, 4, 3), np.array([float(v) for v in lines[-1].split()[:3]])*10


def liquid_path(prep):
    return OUT/f"F_liquid_p{prep}.npz"


def prepare_liquid(prep):
    """Tiled seed-test liquid, remelted with a new velocity seed and re-equilibrated at 220 K (NVT)."""
    import openmm as mm
    from openmm import app, unit
    import system_ice

    if liquid_path(prep).exists():
        return
    xyz, box = read_gro(SEEDTEST/"start.gro")
    shifts = np.array([(i, j, k) for i in range(LIQUID_TILE) for j in range(LIQUID_TILE)
                       for k in range(LIQUID_TILE)])*box
    tiled = np.concatenate([xyz + s for s in shifts])
    gro_path = OUT/f"F_liquid_p{prep}.gro"
    write_gro(gro_path, tiled, box*LIQUID_TILE, f"tiled seed-test liquid, prep {prep}")
    gro = app.GromacsGroFile(str(gro_path))
    top = app.GromacsTopFile(str(gro_path.with_suffix(".top")), periodicBoxVectors=gro.getPeriodicBoxVectors(),
                             includeDir=str(SEEDTEST))
    seed = 7000 + prep
    platform = mm.Platform.getPlatformByName("CUDA")
    system, integrator = system_ice.make_system(top, temp=300.)
    integrator.setRandomNumberSeed(seed)
    ctx = mm.Context(system, integrator, platform, {"Precision": "mixed"})
    ctx.setPeriodicBoxVectors(*gro.getPeriodicBoxVectors())
    ctx.setPositions(gro.getPositions())
    ctx.computeVirtualSites()
    mm.LocalEnergyMinimizer.minimize(ctx, 10., 2000)
    ctx.setVelocitiesToTemperature(300*unit.kelvin, seed)
    integrator.step(LIQUID_MELT_PS*500)
    integrator.setTemperature(220*unit.kelvin)
    integrator.step(LIQUID_EQUIL_PS*500)
    x = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    np.savez_compressed(liquid_path(prep), x=x.reshape(-1, 4, 3), box=box*LIQUID_TILE, velocity_seed=seed)
    gro_path.unlink()
    gro_path.with_suffix(".top").unlink()
    print(f"F_liquid_p{prep}: {len(x)//4} molecules, box {box[0]*LIQUID_TILE:.1f} A", flush=True)


def unwrap_molecules(xyz, box):
    """Keep each molecule whole around its oxygen."""
    return xyz[:, :1] + md_replicas.minimum_image(xyz - xyz[:, :1], box)


def prepare_seed(kind, size, prep, prefix="F"):
    name = f"{prefix}_seed_{kind}_{size}_p{prep}"
    if (OUT/f"{name}.gro").exists():
        return
    rng = np.random.default_rng([2027, KINDS.index(kind), size, prep])
    xyz, box, stack, _, pairs, sequence, dipole = crystal(kind, CRYSTAL, rng)
    xyz = unwrap_molecules(xyz, box)
    o = xyz[:, 0]
    centre = o[rng.integers(len(o))]  # random lattice site; the sphere wraps periodically
    rel = md_replicas.minimum_image(o - centre, box)
    r = np.linalg.norm(rel, axis=1)
    order = np.argsort(r, kind="stable")[:size]
    R = float(r[order].max())
    rotation = Rotation.random(random_state=rng.integers(2**31))
    seed_xyz = rotation.apply((xyz[order] - o[order][:, None] + rel[order][:, None]).reshape(-1, 3)).reshape(-1, 4, 3)
    liquid = np.load(liquid_path(prep))
    lxyz, lbox = liquid["x"], liquid["box"]
    middle = lbox/2
    keep = np.linalg.norm(md_replicas.minimum_image(lxyz[:, 0] - middle, lbox), axis=1) > R + INSERT_MARGIN
    total = np.concatenate([seed_xyz + middle, lxyz[keep]])
    inside = np.zeros(len(o), bool)
    inside[order] = True
    neighbors = [[] for _ in o]
    for a, b in pairs:
        neighbors[a].append(b)
        neighbors[b].append(a)
    interior = np.array([all(inside[j] for j in neighbors[i]) for i in order])
    n = len(total)
    np.savez(OUT/f"{name}.truth.npz", seed=np.r_[np.ones(size, bool), np.zeros(n-size, bool)],
             interior=np.r_[interior, np.zeros(n-size, bool)], stack=np.r_[stack[order], np.zeros(n-size, np.int8)],
             depth=np.r_[R - r[order], np.full(n-size, np.nan)], R=R, centre=middle, restrained=size)
    write_gro(OUT/f"{name}.gro", total, lbox, f"{name} R={R:.1f}A {sequence}")
    meta = dict(kind=kind, size=size, prep=prep, R=R, stacking=sequence, dipole=dipole, molecules=n,
                liquid_removed=int((~keep).sum()), rotation=rotation.as_quat().tolist())
    (OUT/f"{name}.prep.json").write_text(json.dumps(meta, indent=1))


def prepare_slab(kind, prep, name=None, facet="basal"):
    """Periodic slab spanning x and y, liquid above and below in z: two basal (or prism) interfaces.

    For facet="prism" the crystal is built with the stacking axis along z and rotated
    (cyclic x, y, z -> y, z, x) so that the old y axis, normal to {10-10}, becomes the slab normal.
    """
    name = name or f"F_slab_{kind}_p{prep}"
    if (OUT/f"{name}.gro").exists():
        return
    rng = np.random.default_rng([2028 if facet == "basal" else 2029, KINDS.index(kind), prep])
    xyz, box, stack, layer, _, sequence, dipole = crystal(kind, SLAB if facet == "basal" else PRISM, rng)
    xyz = unwrap_molecules(xyz, box)
    if facet == "prism":
        xyz, box = xyz[..., [2, 0, 1]], box[[2, 0, 1]]
    thickness = box[2]
    liquid = np.load(liquid_path(prep))
    lxyz, lbox = liquid["x"], liquid["box"]
    lxyz = lxyz - np.floor(lxyz[:, :1]/lbox)*lbox  # whole molecules, oxygens wrapped into the box
    gap = 42.  # A of liquid between periodic images of the slab
    # liquid cut: oxygens inside the slab's x-y cell and a z window of the gap width
    lo = lxyz[:, 0]
    window = (lo[:, 0] < box[0]) & (lo[:, 1] < box[1]) & (lo[:, 2] < gap - 2*INSERT_MARGIN)
    lsel = lxyz[window] + np.array([0., 0., thickness + INSERT_MARGIN])
    total_box = np.array([box[0], box[1], thickness + gap])
    # The cut liquid is not periodic at the slab's x-y period: drop liquid molecules whose oxygen is
    # closer than OVERLAP to any oxygen already kept (crystal first), across the new periodic boundaries.
    lo = np.mod(lsel[:, 0], total_box)
    near_crystal = cKDTree(np.mod(xyz[:, 0], total_box), boxsize=total_box).query_ball_point(lo, OVERLAP)
    keep = np.array([len(hits) == 0 for hits in near_crystal])
    for i, j in cKDTree(lo, boxsize=total_box).query_pairs(OVERLAP, output_type="ndarray"):
        if keep[i] and keep[j]:
            keep[j] = False
    lsel = lsel[keep]
    total = np.concatenate([xyz, lsel])
    # depth below the nearer basal face, measured from the outermost oxygens as for seeds (R - r)
    z = xyz[:, 0, 2]
    bottom, top = z.min(), z.max()
    depth = np.minimum(z - bottom, top - z)
    n, size = len(total), len(xyz)
    np.savez(OUT/f"{name}.truth.npz", seed=np.r_[np.ones(size, bool), np.zeros(n-size, bool)],
             interior=np.r_[np.ones(size, bool), np.zeros(n-size, bool)],
             stack=np.r_[stack, np.zeros(n-size, np.int8)], depth=np.r_[depth, np.full(n-size, np.nan)],
             R=np.inf, restrained=size)
    write_gro(OUT/f"{name}.gro", total, total_box, f"{name} slab {sequence}")
    meta = dict(kind=kind, prep=prep, facet=facet, stacking=sequence, dipole=dipole, molecules=n, crystal=size,
                liquid=len(lsel), box=total_box.tolist())
    (OUT/f"{name}.prep.json").write_text(json.dumps(meta, indent=1))


def simulate(name):
    import openmm as mm
    from openmm import app, unit
    import system_ice

    gro = app.GromacsGroFile(str(OUT/f"{name}.gro"))
    top = app.GromacsTopFile(str(OUT/f"{name}.top"), periodicBoxVectors=gro.getPeriodicBoxVectors(),
                             includeDir=str(SEEDTEST))
    truth = np.load(OUT/f"{name}.truth.npz")
    seed = 8000 + sum(map(ord, name))
    system, integrator = system_ice.make_system(top)
    restraint = mm.CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
    restraint.addGlobalParameter("k", RESTRAINT)
    for p in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(p)
    start_positions = gro.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    for m in range(int(truth["restrained"])):
        restraint.addParticle(4*m, start_positions[4*m].tolist())
    system.addForce(restraint)
    integrator.setRandomNumberSeed(seed)
    ctx = mm.Context(system, integrator, mm.Platform.getPlatformByName("CUDA"), {"Precision": "mixed"})
    ctx.setPeriodicBoxVectors(*gro.getPeriodicBoxVectors())
    ctx.setPositions(gro.getPositions())
    ctx.computeVirtualSites()
    start = time.time()
    mm.LocalEnergyMinimizer.minimize(ctx, 10., 2000)
    ctx.setVelocitiesToTemperature(220*unit.kelvin, seed)
    integrator.step(RESTRAINED_PS*500)
    ctx.setParameter("k", 0.)
    times = LONG_TIMES if name.startswith("G_long") else tuple(STRIDE_PS*(k+1) for k in range(FRAMES))
    frames, now = [], 0
    for t in times:
        integrator.step((t - now)*500)
        now = t
        frames.append(ctx.getState(getPositions=True).getPositions(asNumpy=True)
                      .value_in_unit(unit.angstrom).astype(np.float32))
    box = np.diag(ctx.getState().getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom))
    np.savez_compressed(OUT/f"{name}.npz", x=np.array(frames), box=box, time_ps=np.array(times))
    print(f"{name}: {len(frames[0])//4} molecules, {time.time()-start:.0f} s", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("systems", nargs="*")
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args()
    if args.list:
        print("\n".join(systems() + confirmation_systems()))
        return
    OUT.mkdir(parents=True, exist_ok=True)
    if args.systems == ["prepare-confirmation"]:
        prep = CONFIRM_PREP
        prepare_liquid(prep)
        for kind in KINDS:
            for size in CONFIRM_SIZES:
                prepare_seed(kind, size, prep, prefix="G")
            prepare_slab(kind, prep, name=f"G_prism_{kind}_p{prep}", facet="prism")
            prepare_slab(kind, prep, name=f"G_long_{kind}_p{prep}")
        return
    if args.systems == ["prepare"]:
        for prep in PREPS:
            prepare_liquid(prep)
            for kind in KINDS:
                prepare_slab(kind, prep)
                for size in SIZES:
                    prepare_seed(kind, size, prep)
        return
    selected = {"all": systems(), "confirmation": confirmation_systems()}.get(" ".join(args.systems), args.systems)
    for name in selected:
        if (OUT/f"{name}.npz").exists():
            print(f"{name}: exists, skipped")
            continue
        simulate(name)


if __name__ == "__main__":
    main()
