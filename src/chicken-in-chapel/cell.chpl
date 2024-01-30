// cell.chpl
// Finite-volume Cell functions for Chicken-in-Chapel.
//
// PJ 2024-01-30

module Cell {
  use Math;
  use IO;
  use Geom;
  use Gas;
  use Flow;
  use Face;

  enum FaceId { iminus, iplus, jminus, jplus, kminus, kplus };

  var FaceIdName = [
    FaceId.iminus => "iminus",
    FaceId.iplus => "iplus",
    FaceId.jminus => "jminus",
    FaceId.jplus => "jplus",
    FaceId.kminus => "kminus",
    FaceId.kplus => "kplus"
  ];

  proc FaceId_from_name(name: string): FaceId throws {
    select name {
      when "iminus" do return FaceId.iminus;
      when "iplus" do return FaceId.iplus;
      when "jminus" do return FaceId.jminus;
      when "jplus" do return FaceId.jplus;
      when "kminus" do return FaceId.kminus;
      when "kplus" do return FaceId.kplus;
      otherwise {
        writeln("Invalid face name: ", name);
        throw new owned Error();
      }
    }
  }

  enum SourceTerms { nosource, manufactured_solution };

  var SourceTermsName = [
    SourceTerms.nosource => "nosource",
    SourceTerms.manufactured_solution => "manufactured_solution"
  ];

  proc SourceTerms_from_name(name: string): SourceTerms throws {
    select name {
      when "nosource" do return SourceTerms.nosource;
      when "none" do return SourceTerms.nosource;
      when "manufactured_solution" do return SourceTerms.manufactured_solution;
      otherwise {
        writeln("Invalid source terms name: ", name);
        throw new owned Error();
      }
    }
  }

  // The IO arrangement is to have the list of items defined by the IOName array,
  // followed by the species mass fractions (#speciesCap of these).
  var IOName = [
    "posx", "posy", "posz", "vol", "p", "T", "rho", "e", "a", "velx", "vely", "velz"
  ];
  // 0       1       2       3      4    5    6      7    8    9       10      11

  proc speciesIOName(i: int): string {
    return try!"massf[%i]".format(i);
  }

  record FVCell {
    var pos: Vector3; // centroid location
    var volume: real;
    var iLength, jLength, kLength;  // for use in interpolation and CFL condition
    var fs: FlowState;
    // We will keep connections to the pieces compising the cell as indices
    // into the block's arrays.
    // Although we probably don't need build and keep this data for the structured grid,
    // it simplifies some of the geometry and update code and may ease the use of
    // unstructured grids at a later date.
    var vtx: [0..#8]int = [0, 0, 0, 0, 0, 0, 0, 0];
    var face: [0..#6]int = [0, 0, 0, 0, 0, 0];

    proc ref IOVarSet(j: int, val: real) throws {
      select j {
        when 0 do pos.x = val;
        when 1 do pos.y = val;
        when 2 do pos.z = val;
        when 3 do volume = val;
        when 4 do fs.gs.p = val;
        when 5 do fs.gs.T = val;
        when 6 do fs.gs.rho = val;
        when 7 do fs.gs.e = val;
        when 8 do fs.gs.a = val;
        when 9 do fs.vel.x = val;
        when 10 do fs.vel.y = val;
        when 11 do fs.vel.z = val;
      }
      const n = IOName.size;
      var i = j - n;
      if i >= speciesCap {
        writeln("Invalid index for IOVar lookup, j=", j);
        throw new owned Error();
      }
      fs.gs.massf[i] = val;
    } // end proc IOVarSet

    proc ref IOVarGet(j: int, val: real): real throws {
      select j {
        when 0 do return pos.x;
        when 1 do return pos.y;
        when 2 do return pos.z;
        when 3 do return volume;
        when 4 do return fs.gs.p;
        when 5 do return fs.gs.T;
        when 6 do return fs.gs.rho;
        when 7 do return fs.gs.e;
        when 8 do return fs.gs.a;
        when 9 do return fs.vel.x;
        when 10 do return fs.vel.y;
        when 11 do return fs.vel.z;
      }
      const n = IOName.size;
      var i = j - n;
      if i >= speciesCap {
        writeln("Invalid index for IOVar lookup, j=", j);
        throw new owned Error();
      }
      return fs.gs.massf[i];
    } // end proc IOVarGet

    proc ref estimate_local_dt(const ref inorm: Vector3,
                               const ref jnorm: Vector3,
                               const ref knorm: Vector3,
                                     cfl: real): real
    {
      // We assume that the cells are (roughly) hexagonal and work with
      // velocities normal to the faces.
      var isignal = iLength/(abs(fs.vel.dot(inorm))+fs.gs.a);
      var jsignal = jLength/(abs(fs.vel.dot(jnorm))+fs.gs.a);
      var ksignal = kLength/(abs(fs.vel.dot(knorm))+fs.gs.a);
      return cfl * min(min(isignal,jsignal),ksignal);
    }

    proc ref addSourceTerms(ref dUdt: [DCon]real, isrc: SourceTerms) {
      select isrc {
        when SourceTerms.nosource do return;
        when SourceTerms.manufactured_solution {
          // [TODO] implement the actual calculation.
          foreach i in DS do dUdt[mass+i] += zero;
          dUdt[xMom] += zero;
          dUdt[yMom] += zero;
          dUdt[zMom] += zero;
          dUdt[totEnergy] += zero;
          return;
        }
      } // end select
    } // end proc addSourceTerms

    // These are the spatial (RHS) terms in the semi-discrete governing equations.
    proc ref eval_dUdt(ref dUdt: [DCon]real,
                       const ref faces: []FVFace, // collection of faces in the block
                       isrc: SourceTerms) {
      const ref fim = faces[face[FaceId.iminus]];
      const ref fip = faces[face[FaceId.iplus]];
      const ref fjm = faces[face[FaceId.jminus]];
      const ref fjp = faces[face[FaceId.jplus]];
      const ref fkm = faces[face[FaceId.kminus]];
      const ref fkp = faces[face[FaceId.kplus]];
      // Introducing local variables for the data helps
      // promote coalesced global memory access on the GPU.
      var area_im = fim.area; const ref F_im = fim.F;
      var area_ip = fip.area; const ref F_ip = fip.F;
      var area_jm = fjm.area; const ref F_jm = fjm.F;
      var area_jp = fjp.area; const ref F_jp = fjp.F;
      var area_km = fkm.area; const ref F_km = fkm.F;
      var area_kp = fkp.area; const ref F_kp = fkp.F;
      //
      var vol_inv = one/volume;
      forall i in DCon {
        // Integrate the fluxes across the interfaces that bound the cell.
        var surface_integral = area_im*F_im[i] - area_ip*F_ip[i]
            + area_jm*F_jm[i] - area_jp*F_jp[i]
            + area_km*F_km[i] - area_kp*F_kp[i];
        // Then evaluate the derivatives of conserved quantity.
        // Note that conserved quantities are stored per-unit-volume.
        dUdt[i] = vol_inv*surface_integral;
      }
      //
      if isrc != SourceTerms.nosource then addSourceTerms(dUdt, isrc);
      return;
    } // end eval_dUdt()

  } // end record FVCell
} // end module Cell
