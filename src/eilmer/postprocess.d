/** postprocess.d
 * Eilmer4 compressible-flow simulation code, postprocessing functions.
 *
 * The role of the post-processing functions is just to pick up data
 * from a previously-run simulation and either write plotting files
 * or extract interesting pieces of data.
 *
 * Author: Peter J. and Rowan G.
 * First code: 2015-06-09
 */

module postprocess;

import std.math;
import std.stdio;
import std.conv;
import std.format;
import std.string;
import std.regex;
import std.algorithm;
import std.bitmanip;
import std.stdint;
import std.range;
import nm.complex;
import nm.number;
import gzip;
import fileutil;
import geom;
import gas;
import globalconfig;
import flowsolution;
import solidsolution;
import vtk_writer;
import tecplot_writer;
import tecplot_writer_classic;


void post_process(string plotDir, bool listInfoFlag, string tindxPlot,
                  string addVarsStr, string luaRefSoln,
                  bool vtkxmlFlag, bool binary_format,
                  bool tecplotBinaryFlag, bool tecplotAsciiFlag, bool tecplotAsciiLegacyFlag,
                  string outputFileName, string sliceListStr,
                  string surfaceListStr, string extractStreamStr, string trackWaveStr,
                  string extractLineStr, string computeLoadsOnGroupStr,
                  string probeStr, string outputFormat,
                  string normsStr, string regionStr,
                  string extractSolidLineStr, string tag)
{
    read_config_file();
    string jobName = GlobalConfig.base_file_name;
    string outName = jobName;
    if (GlobalConfig.new_flow_format){
        outName ~= "." ~ tag;
    }
    //
    string[] addVarsList;
    addVarsStr = addVarsStr.strip();
    addVarsStr = addVarsStr.replaceAll(regex("\""), "");
    if (addVarsStr.length > 0) {
        addVarsList = addVarsStr.split(",");
    }
    //
    auto times_dict = readTimesFile(jobName);
    auto tindx_list = times_dict.keys;
    sort(tindx_list);
    int[] tindx_list_to_plot;
    tindxPlot = tindxPlot.replaceAll(regex("\""), "");
    switch (tindxPlot) {
    case "all":
        tindx_list_to_plot = tindx_list.dup;
        break;
    case "9999":
    case "last":
        tindx_list_to_plot ~= tindx_list[$-1];
        break;
    default:
        // We assume that the command-line argument was an integer
        // or a comma-separated list of integers.
        auto items = tindxPlot.split(",");
        foreach (item; items) {
            if (item.length > 0) {
                int i = to!int(item);
                if ((i >= tindx_list[0]) && (i <= tindx_list[$-1])) {
                    tindx_list_to_plot ~= i;
                }
            }
        }
        sort(tindx_list_to_plot);
    } // end switch
    //
    if (listInfoFlag) {
        writeln("Some information about this simulation.");
        writeln("  nFluidBlocks= ", GlobalConfig.nFluidBlocks);
        writeln("  nSolidBlocks= ", GlobalConfig.nSolidBlocks);
        writeln("  last tindx= ", tindx_list[$-1]);
        writeln("  Flow Variables:");
        // Dip into the top of a solution file that is likely to be present
        // to get the variable names, as saved by the simulation.
        double sim_time;
        switch (GlobalConfig.flow_format) {
        case "gziptext": goto default;
        case "rawbinary": throw new Error("not yet implemented PJ 2017-09-02");
        default:
            string fileName = make_file_name!"flow"(jobName, to!int(0), 0, "gz");
            auto byLine = new GzipByLine(fileName);
            auto line = byLine.front; byLine.popFront();
            formattedRead(line, " %g", &sim_time);
            line = byLine.front; byLine.popFront();
            auto variableNames = line.strip().split();
            foreach (ref var; variableNames) { var = var.replaceAll(regex("\""), ""); }
            foreach (i; 0 .. variableNames.length) {
                writeln(format("%4d %s", i, variableNames[i]));
            }
        } // end switch flow_format
        if ( GlobalConfig.nSolidBlocks > 0 ) {
            writeln("  Solid Variables:");
            // Dip into the top of a solid solution file that is
            // likely to be present to get the variable names
            // as saved by the simulation.
            string fileName = make_file_name!"solid"(jobName, 0, 0, "gz");
            auto byLine = new GzipByLine(fileName);
            auto line = byLine.front; byLine.popFront();
            formattedRead(line, " %g", &sim_time);
            line = byLine.front; byLine.popFront();
            auto variableNames = line.strip().split();
            foreach (ref var; variableNames) { var = var.replaceAll(regex("\""), ""); }
            foreach (i; 0 .. variableNames.length) {
                writeln(format("%4d %s", i, variableNames[i]));
            }
        } // end if nSolidBlocks > 0
    } // end if listInfoFlag
    //
    if (vtkxmlFlag) {
        ensure_directory_is_present(plotDir);
        //
        writeln("writing flow-solution VTK-XML files to directory \"", plotDir, "\"");
        File visitFile = File(plotDir~"/"~outName~".visit", "w");
        // For each time index, the visit justs lists the names of the files for individual blocks.
        visitFile.writef("!NBLOCKS %d\n", GlobalConfig.nFluidBlocks);
        File pvdFile = begin_PVD_file(plotDir~"/"~outName~".pvd");
        foreach (tindx; tindx_list_to_plot) {
            writeln("  tindx= ", tindx);
            auto soln = new FlowSolution(jobName, ".", tindx, GlobalConfig.nFluidBlocks, -1, GlobalConfig.flow_format, tag);
            soln.add_aux_variables(addVarsList, tag);
            if (luaRefSoln.length > 0) soln.subtract_ref_soln(luaRefSoln);
            string pvtuFileName = outName~format("-t%04d", tindx)~".pvtu";
            File pvtuFile = begin_PVTU_file(plotDir~"/"~pvtuFileName, soln.flowBlocks[0].variableNames);
            foreach (jb; 0 .. soln.nBlocks) {
                string vtuFileName = outName~format("-b%04d-t%04d.vtu", jb, tindx);
                add_dataset_to_PVD_file(pvdFile, soln.sim_time, vtuFileName);
                add_piece_to_PVTU_file(pvtuFile, vtuFileName);
                visitFile.writef("%s\n", vtuFileName);
                write_VTU_file(soln.flowBlocks[jb], soln.gridBlocks[jb], plotDir~"/"~vtuFileName, binary_format);
            }
            finish_PVTU_file(pvtuFile);
        } // foreach tindx
        finish_PVD_file(pvdFile);
        visitFile.close();
        //
        if ( GlobalConfig.nSolidBlocks > 0 ) {
            writeln("writing solid VTK-XML files to directory \"", plotDir, "\"");
            visitFile = File(plotDir~"/"~jobName~"-solid.visit", "w");
            // For each time index, the visit justs lists the names of the files for individual blocks.
            visitFile.writef("!NBLOCKS %d\n", GlobalConfig.nSolidBlocks);
            pvdFile = begin_PVD_file(plotDir~"/"~jobName~"-solid.pvd");
            foreach (tindx; tindx_list_to_plot) {
                writeln("  tindx= ", tindx);
                auto soln = new SolidSolution(jobName, ".", tindx, GlobalConfig.nSolidBlocks);
                if (luaRefSoln.length > 0) soln.subtract_ref_soln(luaRefSoln);
                string pvtuFileName = jobName~format("-solid-t%04d", tindx)~".pvtu";
                File pvtuFile = begin_PVTU_file(plotDir~"/"~pvtuFileName, soln.solidBlocks[0].variableNames, false);
                foreach (jb; 0 .. soln.nBlocks) {
                    string vtuFileName = jobName~format("-solid-b%04d-t%04d.vtu", jb, tindx);
                    add_piece_to_PVTU_file(pvtuFile, vtuFileName);
                    add_dataset_to_PVD_file(pvdFile, soln.sim_time, vtuFileName);
                    visitFile.writef("%s\n", vtuFileName);
                    write_VTU_file(soln.solidBlocks[jb], soln.gridBlocks[jb], plotDir~"/"~vtuFileName, binary_format);
                }
                finish_PVTU_file(pvtuFile);
            } // foreach tindx
            finish_PVD_file(pvdFile);
            visitFile.close();
        } // end if nSolidBlocks > 0
    } // end if vtkxml
    //
    version(with_tecplot_binary) {
    if (tecplotBinaryFlag) {
        // Use Pierpaolo's Tecplot-writer module.
        ensure_directory_is_present(plotDir);
        writeln("Writing Tecplot (binary) file(s) to directory \"", plotDir, "\"");
        foreach (tindx; tindx_list_to_plot) {
            writeln("  tindx= ", tindx);
            double timeStamp = times_dict[tindx];
            auto soln = new FlowSolution(jobName, ".", tindx, GlobalConfig.nFluidBlocks, -1, GlobalConfig.flow_format, tag);
            soln.add_aux_variables(addVarsList, tag);
            string fname = format("%s/%s.t%04d", plotDir, jobName, tindx);
            if ( writeTecplotBinaryHeader(jobName, tindx, fname, soln.flowBlocks[0].variableNames) != 0 ) {
                string errMsg = format("Tecplot binary output failed for tindx: %d", tindx);
                throw new FlowSolverException(errMsg);
            }
            foreach (jb; 0 .. GlobalConfig.nFluidBlocks) {
                int zoneType;
                size_t[][] connList;
                prepareGridConnectivity(soln.gridBlocks[jb], zoneType, connList);
                if ( writeTecplotBinaryZoneHeader(soln.flowBlocks[jb], soln.gridBlocks[jb], jb,
                                                  soln.flowBlocks[jb].variableNames, timeStamp, zoneType) != 0 ) {
                    string errMsg = format("Tecplot binary output failed for block: %d", jb);
                    throw new FlowSolverException(errMsg);
                }
                writeTecplotBinaryZoneData(soln.flowBlocks[jb], soln.gridBlocks[jb],
                                           soln.flowBlocks[jb].variableNames, connList);
            }
            if ( closeTecplotBinaryFile() != 0 ) {
                string errMsg = format("Closing of Tecplot binary file failed for tindx: %d", tindx);
                throw new FlowSolverException(errMsg);
            }
        }
    } // end if tecplotBinaryFlag
    }
    else {
    if (tecplotBinaryFlag) {
        string errMsg = "This version of e4shared was NOT compiled with support for Tecplot binary files.";
        throw new FlowSolverException(errMsg);
    }
    }
    if (tecplotAsciiFlag) {
        // Use Pierpaolo's tecplot-writer module.
        ensure_directory_is_present(plotDir);
        writeln("writing Tecplot ASCII-text file(s) to directory \"", plotDir, "\"");
        foreach (tindx; tindx_list_to_plot) {
            writeln("  tindx= ", tindx);
            double timeStamp = times_dict[tindx];
            auto soln = new FlowSolution(jobName, ".", tindx, GlobalConfig.nFluidBlocks, -1, GlobalConfig.flow_format, tag);
            soln.add_aux_variables(addVarsList, tag);
            string fname = format("%s/%s-t%04d.tec", plotDir, jobName, tindx);
            File fp = writeTecplotAsciiHeader(jobName, tindx, fname, soln.flowBlocks[0].variableNames);
            foreach (jb; 0 .. GlobalConfig.nFluidBlocks) {
                int zoneType;
                size_t[][] connList;
                prepareGridConnectivity(soln.gridBlocks[jb], zoneType, connList);
                writeTecplotAsciiZoneHeader(soln.flowBlocks[jb], soln.gridBlocks[jb], jb, fp,
                                            soln.flowBlocks[jb].variableNames, timeStamp, zoneType);
                writeTecplotAsciiZoneData(soln.flowBlocks[jb], soln.gridBlocks[jb], fp,
                                          soln.flowBlocks[jb].variableNames, connList);
            }
            fp.close();
        }
    } // end
    if (tecplotAsciiLegacyFlag) {
        // For structured-grid blocks, we had an old Tecplot writer from long ago.
        ensure_directory_is_present(plotDir);
        writeln("writing Tecplot ASCII-text (Legacy) file(s) to directory \"", plotDir, "\"");
        foreach (tindx; tindx_list_to_plot) {
            writeln("  tindx= ", tindx);
            auto soln = new FlowSolution(jobName, ".", tindx, GlobalConfig.nFluidBlocks, -1, GlobalConfig.flow_format, tag);
            // Temporary check for unstructured grids. Remove when unstructured version is implemented.
            foreach ( grid; soln.gridBlocks ) {
                if ( grid.grid_type == Grid_t.unstructured_grid ) {
                    throw new FlowSolverException("Tecplot output not currently available for unstructured grids.");
                }
            }
            soln.add_aux_variables(addVarsList, tag);
            if (luaRefSoln.length > 0) soln.subtract_ref_soln(luaRefSoln);
            auto t = times_dict[tindx];
            write_Tecplot_file(jobName, plotDir, soln, tindx);
        } // foreach tindx
        if ( GlobalConfig.nSolidBlocks > 0 ) {
            throw new FlowSolverException("Tecplot output not currently available for solid blocks.");
        //     writeln("writing solid Tecplot file(s) to directory \"", plotDir, "\"");
        //     foreach (tindx; tindx_list_to_plot) {
        //      writeln("  tindx= ", tindx);
        //      auto soln = new SolidSolution(jobName, ".", tindx, GlobalConfig.nSolidBlocks);
        //      if (luaRefSoln.length > 0) soln.subtract_ref_soln(luaRefSoln);
        //      write_Tecplot_file(jobName, plotDir, soln, tindx);
        //     } // foreach tindx
        } // end if nSolidBlocks > 0
    } // end if tecplot
    //
    if (probeStr.length > 0) {
        writeln("Probing flow solution at specified points.");
        // The output may go to a user-specified file, or stdout.
        File outFile;
        if (outputFileName.length > 0) {
            outFile = File(outputFileName, "w");
            writeln("Output will be sent to File: ", outputFileName);
        } else {
            outFile = stdout;
        }
        probeStr = probeStr.strip();
        probeStr = probeStr.replaceAll(regex("\""), "");
        double[] xp, yp, zp;
        foreach(triple; probeStr.split(";")) {
            auto items = triple.split(",");
            xp ~= to!double(items[0]);
            yp ~= to!double(items[1]);
            zp ~= to!double(items[2]);
        }
        foreach (tindx; tindx_list_to_plot) {
            writeln("  tindx= ", tindx);
            auto soln = new FlowSolution(jobName, ".", tindx, GlobalConfig.nFluidBlocks, -1, GlobalConfig.flow_format, tag);
            soln.add_aux_variables(addVarsList, tag);
            if (luaRefSoln.length > 0) soln.subtract_ref_soln(luaRefSoln);
            if (outputFormat == "gnuplot") {
                outFile.writeln(soln.flowBlocks[0].variable_names_as_string(true));
            }
            foreach (ip; 0 .. xp.length) {
                auto nearest = soln.find_nearest_cell_centre(xp[ip], yp[ip], zp[ip]);
                size_t ib = nearest[0]; size_t i = nearest[1];
                if (outputFormat == "gnuplot") {
                    outFile.writeln(soln.flowBlocks[ib].values_as_string(i));
                } else {
                    // Assume that pretty format was requested.
                    outFile.writefln("Block[%d], cell[%d]:", ib, i);
                    outFile.writefln("  pos=(%s, %s, %s)m, volume=%s m^^3",
                                     soln.get_value_str(ib, i, "pos.x"), soln.get_value_str(ib, i, "pos.y"),
                                     soln.get_value_str(ib, i, "pos.z"), soln.get_value_str(ib, i, "volume"));
                    outFile.writefln("  pos=(%s, %s, %s)m, volume=%s m^^3",
                                     soln.get_value_str(ib, i, "pos.x"), soln.get_value_str(ib, i, "pos.y"),
                                     soln.get_value_str(ib, i, "pos.z"), soln.get_value_str(ib, i, "volume"));
                    outFile.writefln("  rho=%s kg/m^^3, p=%s Pa, T=%s K, u=%s J/kg",
                                     soln.get_value_str(ib, i, "rho"), soln.get_value_str(ib, i, "p"),
                                     soln.get_value_str(ib, i, "T"), soln.get_value_str(ib, i, "u"));
                    outFile.writefln("  vel=(%s, %s, %s)m/s, a=%s m/s",
                                     soln.get_value_str(ib, i, "vel.x"), soln.get_value_str(ib, i, "vel.y"),
                                     soln.get_value_str(ib, i, "vel.z"), soln.get_value_str(ib, i, "a"));
                    outFile.writefln("  M_local=%s, pitot_p=%s Pa, total_p=%s Pa, total_h=%s J/kg",
                                     soln.get_value_str(ib, i, "M_local"), soln.get_value_str(ib, i, "pitot_p"),
                                     soln.get_value_str(ib, i, "total_p"), soln.get_value_str(ib, i, "total_h"));
                    outFile.writefln("  mu=%s Pa.s, k=%s W/(m.K)", soln.get_value_str(ib, i, "mu"),
                                     soln.get_value_str(ib, i, "k"));
                    outFile.writefln("  mu_t=%s Pa.s, k_t=%s W/(m.K), tke=%s (m/s)^^2, omega=%s 1/s",
                                     soln.get_value_str(ib, i, "mu_t"), soln.get_value_str(ib, i, "k_t"),
                                     soln.get_value_str(ib, i, "tke"), soln.get_value_str(ib, i, "omega"));
                    outFile.writefln("  massf=[%s]", soln.get_massf_str(ib, i));
                }
            }
        } // end foreach tindx

        if (GlobalConfig.nSolidBlocks > 0) {
            foreach (tindx; tindx_list_to_plot) {
                writeln("  tindx= ", tindx);
                auto soln = new SolidSolution(jobName, ".", tindx, GlobalConfig.nSolidBlocks);
                if (luaRefSoln.length > 0) soln.subtract_ref_soln(luaRefSoln);
                if (outputFormat == "gnuplot") {
                    outFile.writeln(soln.solidBlocks[0].variable_names_as_string());
                }
                foreach (ip; 0 .. xp.length) {
                    auto nearest = soln.find_nearest_cell_centre(xp[ip], yp[ip], zp[ip]);
                    size_t ib = nearest[0]; size_t i = nearest[1]; size_t j = nearest[2]; size_t k = nearest[3];
                    if (outputFormat == "gnuplot") {
                        outFile.writeln(soln.solidBlocks[ib].values_as_string(i, j, k));
                    } else {
                        // Assume that pretty format was requested.
                        outFile.writefln("SolidBlock[%d], cell[%d]:", ib, i);
                        outFile.writefln("  pos=(%s, %s, %s)m, volume=%s m^^3",
                                         soln.get_value_str(ib, i, j, k, "pos.x"), soln.get_value_str(ib, i, j, k, "pos.y"),
                                         soln.get_value_str(ib, i, j, k, "pos.z"), soln.get_value_str(ib, i, j, k, "volume"));
                        outFile.writefln("  e=%s J/kg, T=%s K",
                                         soln.get_value_str(ib, i, j, k, "e"), soln.get_value_str(ib, i, j, k, "T"));
                    }
                }
            } // end foreach tindx
        } // end if nSolidBlocks > 0
    } // end if probeStr
    //
    if (sliceListStr.length > 0) {
        writeln("Extracting slices of the flow solution.");
        // The output may go to a user-specified file, or stdout.
        File outFile;
        if (outputFileName.length > 0) {
            outFile = File(outputFileName, "w");
            writeln("Output will be sent to File: ", outputFileName);
        } else {
            outFile = stdout;
        }
        bool header_written = false;
        foreach (tindx; tindx_list_to_plot) {
            writeln("  tindx= ", tindx);
            auto soln = new FlowSolution(jobName, ".", tindx, GlobalConfig.nFluidBlocks, -1, GlobalConfig.flow_format, tag);
            soln.add_aux_variables(addVarsList, tag);
            if (luaRefSoln.length > 0) soln.subtract_ref_soln(luaRefSoln);
            //
            if (!header_written) {
                // Gnuplot column labels.
                outFile.writeln(soln.flowBlocks[0].variable_names_as_string(true, false, true));
                header_written = true;
            } else {
                // Gnuplot datasets are separated by two blank lines.
                outFile.write("\n\n");
            }
            foreach (sliceStr; sliceListStr.split(";")) {
                auto rangeStrings = sliceStr.split(",");
                auto blk_range = decode_range_indices(rangeStrings[0], 0, soln.nBlocks);
                foreach (ib; blk_range[0] .. blk_range[1]) {
                    auto blk = soln.flowBlocks[ib];
                    // We need to do the decode in the context of each block because
                    // the upper limits to the indices are specific to the block.
                    auto i_range = decode_range_indices(rangeStrings[1], 0, blk.nic);
                    auto j_range = decode_range_indices(rangeStrings[2], 0, blk.njc);
                    auto k_range = decode_range_indices(rangeStrings[3], 0, blk.nkc);
                    foreach (k; k_range[0] .. k_range[1]) {
                        foreach (j; j_range[0] .. j_range[1]) {
                            foreach (i; i_range[0] .. i_range[1]) {
                                outFile.writeln(blk.values_as_string(i,j,k));
                            }
                        }
                    }
                } // end foreach ib
            } // end foreach sliceStr
        } // end foreach tindx
    } // end if sliceListStr
    //
    if (surfaceListStr.length > 0) {
        writeln("Extracting named surfaces of the flow solution.");
        writeln("writing VTK-XML files to directory \"", plotDir, "\"");
        ensure_directory_is_present(plotDir);
        string surfaceCollectionName = outputFileName;
        if (surfaceCollectionName.length == 0) {
            throw new Exception("Expected name for surface collection to be provided with --output-file");
        }
        File pvdFile = begin_PVD_file(plotDir~"/"~surfaceCollectionName~".pvd");
        foreach (tindx; tindx_list_to_plot) {
            writeln("  tindx= ", tindx);
            auto soln = new FlowSolution(jobName, ".", tindx, GlobalConfig.nFluidBlocks, -1,
                                         GlobalConfig.flow_format, tag);
            soln.add_aux_variables(addVarsList, tag);
            string pvtuFileName = surfaceCollectionName~format("-t%04d", tindx)~".pvtu";
            File pvtuFile = begin_PVTU_file(plotDir~"/"~pvtuFileName, soln.flowBlocks[0].variableNames);
            foreach (surfaceStr; surfaceListStr.split(";")) {
                auto itemStrings = surfaceStr.split(",");
                size_t blk_indx = to!size_t(itemStrings[0]);
                string boundary_id = itemStrings[1];
                size_t boundary_indx;
                if (canFind(face_name, boundary_id)) {
                    boundary_indx = face_index(boundary_id);
                } else {
                    boundary_indx = to!size_t(boundary_id);
                }
                writefln("Assemble boundary surface %d for block %d", boundary_indx, blk_indx);
                string vtuFileName = format("%s-t%04d-blk-%04d-surface-%s.vtu",
                                            surfaceCollectionName, tindx, blk_indx, boundary_id);
                auto surf_grid = soln.gridBlocks[blk_indx].get_boundary_grid(boundary_indx);
                size_t[] surf_cells = soln.gridBlocks[blk_indx].get_list_of_boundary_cells(boundary_indx);
                if (surf_cells.length > 0) {
                    size_t new_dimensions = surf_grid.dimensions;
                    size_t new_nic;
                    size_t new_njc;
                    size_t new_nkc;
                    final switch (surf_grid.grid_type) {
                    case Grid_t.structured_grid:
                        new_nic = max(surf_grid.niv-1, 1);
                        new_njc = max(surf_grid.njv-1, 1);
                        new_nkc = max(surf_grid.nkv-1, 1);
                        if (new_nic*new_njc*new_nkc != surf_cells.length) {
                            throw new Exception("Mismatch in number of cells for new surface grid.");
                        }
                        break;
                    case Grid_t.unstructured_grid:
                        new_nic = surf_cells.length;
                        new_njc = 1;
                        new_nkc = 1;
                        break;
                    }
                    auto surf_flow = new FluidBlockLite(soln.flowBlocks[blk_indx], surf_cells,
                                                        new_dimensions, new_nic, new_njc, new_nkc);
                    add_dataset_to_PVD_file(pvdFile, soln.sim_time, vtuFileName);
                    add_piece_to_PVTU_file(pvtuFile, vtuFileName);
                    write_VTU_file(surf_flow, surf_grid, plotDir~"/"~vtuFileName, binary_format);
                } else {
                    writeln("Warning: no flow cells associated with the new surface grid.");
                }
            } // end foreach surfaceStr
            finish_PVTU_file(pvtuFile);
        } // end foreach tindx
        finish_PVD_file(pvdFile);
    } // end if surfaceListStr
    //
    if (extractStreamStr.length > 0) {
        writeln("Extracting data along a streamline of the flow solution.");
        // The output may go to a user-specified file, or stdout.
        File outFile;
        if (outputFileName.length > 0) {
            outFile = File(outputFileName, "w");
            writeln("Output will be sent to File: ", outputFileName);
        } else {
            outFile = stdout;
        }
        extractStreamStr = extractStreamStr.strip();
        extractStreamStr = extractStreamStr.replaceAll(regex("\""), "");
        double[] xp, yp, zp;
        foreach(triple; extractStreamStr.split(";")) {
            auto items = triple.split(",");
            xp ~= to!double(items[0]);
            yp ~= to!double(items[1]);
            zp ~= to!double(items[2]);
        }
        double stepSize = 1e-06; // set a temporal step size
        double xInit, yInit, zInit;
        double xOld, yOld, zOld;
        double xNew, yNew, zNew;
        foreach (tindx; tindx_list_to_plot) {
            writeln("  tindx= ", tindx);
            auto soln = new FlowSolution(jobName, ".", tindx, GlobalConfig.nFluidBlocks, -1,
                                         GlobalConfig.flow_format, tag);
            soln.add_aux_variables(addVarsList, tag);
            if (luaRefSoln.length > 0) soln.subtract_ref_soln(luaRefSoln);
            outFile.writeln("# xStreamPos ", "yStreamPos ", "zStreamPos ", "relDistance ",
                            soln.flowBlocks[0].variable_names_as_string());
            foreach (ip; 0 .. xp.length) {
                outFile.writeln("# streamline locus point: ", xp[ip], ", ", yp[ip], ", ", zp[ip]);
                auto identity = soln.find_enclosing_cell(xp[ip], yp[ip], zp[ip]);
                size_t ib = identity[0]; size_t idx = identity[1]; size_t found = identity[2];
                if (found == 0) { // out of domain bounds
                    writeln("User defined point not in solution domain bounds");
                    break;
                }
                else { // store initial cell data
                    xInit = soln.flowBlocks[ib]["pos.x", idx];
                    yInit = soln.flowBlocks[ib]["pos.y", idx];
                    zInit = soln.flowBlocks[ib]["pos.z", idx];
                    outFile.writeln(xInit, " ", yInit, " ", zInit, " 0 ",
                                    soln.flowBlocks[ib].values_as_string(idx));
                }
                // we need to travel both forward (direction = 1) and backward (direction = -1)
                int[] direction = [-1, 1];
                double min = 1e-6;
                foreach (direct; direction) {
                    found = 1;
                    xOld = xInit; yOld = yInit; zOld = zInit;
                    double distance = 0.0; // relative distance along streamline
                    while (found == 1) { // while we have a cell in the domain
                        double vx = soln.flowBlocks[ib]["vel.x", idx];
                        double vy = soln.flowBlocks[ib]["vel.y", idx];
                        double vz = soln.flowBlocks[ib]["vel.z", idx];
                        double dx = direct*vx*stepSize;
                        double dy = direct*vy*stepSize;
                        double dz = direct*vz*stepSize;
                        distance += direct*sqrt(dx*dx + dy*dy + dz*dz);
                        xNew = xOld + dx; yNew = yOld + dy; zNew = zOld + dz;
                        identity = soln.find_enclosing_cell(xNew, yNew, zNew);
                        if (identity[0] == ib && identity[1] == idx) {
                            // did not step outside current cell
                            stepSize = stepSize*2.0; found = identity[2];
                        } else {
                            ib = identity[0]; idx = identity[1]; found = identity[2];
                            if (found == 1) {
                                outFile.writeln(xNew, " ", yNew, " ", zNew, " ", distance, " ",
                                                soln.flowBlocks[ib].values_as_string(idx));
                            }
                            xOld = xNew; yOld = yNew; zOld = zNew;
                            stepSize = 1e-06;
                        } // end else
                    } // end while
                } // end foreach direction
            } // end for each xp.length
        } // end foreach tindx
    } // end if streamlineStr
    //
    if (trackWaveStr.length > 0) {
        writeln("Tracking a wave through the flow solution.");
        // The output may go to a user-specified file, or stdout.
        File outFile;
        if (outputFileName.length > 0) {
            outFile = File(outputFileName, "w");
            writeln("Output will be sent to File: ", outputFileName);
        } else {
            outFile = stdout;
        }
        trackWaveStr = trackWaveStr.strip();
        trackWaveStr = trackWaveStr.replaceAll(regex("\""), "");
        double[] xp, yp, zp;
        Vector3[] SliceNormal;
        // Extract starting point coordinates and normal vector of the slice
        foreach(pointAndNormal; trackWaveStr.split(";")) {
            auto items = pointAndNormal.split(",");
            xp ~= to!double(items[0]);
            yp ~= to!double(items[1]);
            zp ~= to!double(items[2]);

            if (items.length == 6) {
                SliceNormal ~= Vector3(to!double(items[3]),to!double(items[4]),to!double(items[5]));
            } else {
                // if no normal vector is given assume a plane parallel to x-y
                SliceNormal ~= Vector3(0.0,0.0,1.0);
            }
            SliceNormal[$-1] = unit(SliceNormal[$-1]);
        }
        double stepSize = 1e-06; // set a temporal step size
        double xInit, yInit, zInit;
        foreach (tindx; tindx_list_to_plot) {
            writeln("  tindx= ", tindx);
            auto soln = new FlowSolution(jobName, ".", tindx, GlobalConfig.nFluidBlocks, -1,
                                         GlobalConfig.flow_format, tag);
            soln.add_aux_variables(addVarsList, tag);
            if (luaRefSoln.length > 0) soln.subtract_ref_soln(luaRefSoln);
            outFile.writeln("# xWavePos ", "yWavePos ", "zWavePos ", "relDistance ",
                            soln.flowBlocks[0].variable_names_as_string());
            foreach (ip; 0 .. xp.length) {
                outFile.writeln("# wave locus point: ", xp[ip], ", ", yp[ip], ", ", zp[ip]);
                auto identity = soln.find_enclosing_cell(xp[ip], yp[ip], zp[ip]);
                size_t ib = identity[0]; size_t idx = identity[1]; size_t found = identity[2];
                if (found == 0) { // out of domain bounds
                    writeln("User defined point not in solution domain bounds");
                    break;
                }
                else { // store initial cell data
                    xInit = soln.flowBlocks[ib]["pos.x", idx];
                    yInit = soln.flowBlocks[ib]["pos.y", idx];
                    zInit = soln.flowBlocks[ib]["pos.z", idx];
                    xInit = xp[ip];
                    yInit = yp[ip];
                    zInit = zp[ip];
                    outFile.writeln(xInit, " ", yInit, " ", zInit, " ",
                                    soln.flowBlocks[ib].values_as_string(idx));
                }
                size_t ibInit = ib; size_t idxInit = idx;
                // we need to travel both forward (direction = 1) and backward (direction = -1)
                int[] direction = [-1, 1];
                double min = 1e-6;
                foreach (direct; direction) {
                    // In every slice there are two waves emanating from every one point
                    foreach (dir; direction) {
                        found = 1;
                        Vector3 P0 = Vector3(xInit,yInit,zInit);
                        number distance = 0.0; // relative distance along streamline
                        outFile.writeln("# New Wave");
                        ib = ibInit; idx = idxInit;
                        while (found == 1) { // while we have a cell in the domain
                            Vector3 vlocal = Vector3(soln.flowBlocks[ib]["vel.x", idx],
                                                     soln.flowBlocks[ib]["vel.y", idx],
                                                     soln.flowBlocks[ib]["vel.z", idx]);
                            Vector3 dStream = vlocal*direct*stepSize;
                            Vector3 P1 = P0 + dStream;

                            // define slice as n1*x+n2*y+n3*z+sliceConst = 0
                            number sliceConst = -dot(SliceNormal[ip],P0);
                            number coneSliceDist = abs(dot(SliceNormal[ip],P1)+sliceConst);

                            // calculate local Mach number and angle
                            number alocal = soln.flowBlocks[ib]["a", idx];
                            number Mlocal = geom.abs(vlocal)/alocal;
                            number MachAngle;

                            // check if flow is supersonic
                            if (Mlocal >= 1) {
                                MachAngle = asin(1/Mlocal);
                            } else {
                                writeln("Subsonic flow encountered");
                                break;
                            }

                            Vector3 P2 = P1+(SliceNormal[ip]*coneSliceDist);

                            // check if Point2 is on the specified slice if not flip normal vector
                            if (abs(dot(SliceNormal[ip],P2)+sliceConst) >= 1e-16) {
                                P2 = P1-(SliceNormal[ip]*coneSliceDist);
                            }

                            // calculate angle in between the slice and the velocity vector
                            Vector3 SliceVec = P2-P0;
                            number beta = acos(dot(SliceVec,dStream)/(geom.abs(SliceVec)*geom.abs(dStream)));

                            // Check if the slice intersects the Mach cone
                            if (beta > MachAngle) {
                                writeln("Specified slice doesn't intersect the Mach cone");
                                writeln("Locus Point: ", xp[ip], ", ", yp[ip], ", ", zp[ip]);
                                writeln("Slice Normal: ", SliceNormal[ip].x, ", ",
                                        SliceNormal[ip].y, ", ", SliceNormal[ip].z);
                                writeln("Direction: ", direct);
                                break;
                            }

                            // P2 is closer to P0 than P1, need to calculate the Mach cone
                            // at the new position.
                            // Project vector P0P2 onto the streamline.
                            dStream = unit(dStream);
                            Vector3 SliceVecProj = dot(SliceVec,dStream)*dStream;
                            Vector3 P3 = P0+SliceVecProj;

                            // Calculate radius of the Mach cone at the new point
                            number rCone = tan(MachAngle)*geom.abs(SliceVecProj);
                            Vector3 P2P3 = P3-P2;
                            number dP2P3 = geom.abs(P2P3);
                            number dP2P4 = sqrt(rCone^^2-dP2P3^^2);

                            Vector3 SliceVec2 = cross(SliceVec,SliceNormal[ip]);
                            SliceVec2 = unit(SliceVec2);
                            Vector3 P4 = P2+(dir*dP2P4*SliceVec2);

                            // Calculate direction and length of the wave segment
                            Vector3 Wave = P4-P0;
                            distance += direct*geom.abs(Wave);

                            number WaveAngle = acos(dot(Wave,dStream)/(geom.abs(Wave)*geom.abs(dStream)));

                            identity = soln.find_enclosing_cell(P4.x.re, P4.y.re, P4.z.re);
                            if (identity[0] == ib && identity[1] == idx) {
                                // did not step outside current cell
                                stepSize = stepSize*2.0; found = identity[2];
                                writeln("Increasing Step size");
                            } else {
                                ib = identity[0]; idx = identity[1]; found = identity[2];
                                if (found == 1) {
                                    outFile.writeln(P4.x, " ", P4.y, " ", P4.z, " ", distance, " ",
                                                    soln.flowBlocks[ib].values_as_string(idx));
                                }
                                P0 = P4;
                                stepSize = 1e-06;
                            } // end else
                         } // end while
                    } // end foreach wave direction
                } // end foreach direction
            } // end for each xp.length
        } // end foreach tindx
    } // end if trackWaveStr
    //
    if (extractLineStr.length > 0) {
        writeln("Extracting data along a straight line between end points.");
        // The output may go to a user-specified file, or stdout.
        File outFile;
        if (outputFileName.length > 0) {
            outFile = File(outputFileName, "w");
            writeln("Output will be sent to File: ", outputFileName);
        } else {
            outFile = stdout;
        }
        bool header_written = false;
        extractLineStr = extractLineStr.strip();
        extractLineStr = extractLineStr.replaceAll(regex("\""), "");
        foreach (tindx; tindx_list_to_plot) {
            writeln("  tindx= ", tindx);
            auto soln = new FlowSolution(jobName, ".", tindx, GlobalConfig.nFluidBlocks, -1,
                                         GlobalConfig.flow_format, tag);
            soln.add_aux_variables(addVarsList, tag);
            if (luaRefSoln.length > 0) soln.subtract_ref_soln(luaRefSoln);
            if (!header_written) {
                // Gnuplot column labels.
                outFile.writeln(soln.flowBlocks[0].variable_names_as_string(true));
                header_written = true;
            } else {
                // Gnuplot 2 blank lines separates datasets.
                outFile.write("\n\n");
            }
            size_t[2][] cells_found; // accumulate the identies of the cells found here
            foreach(lineStr; extractLineStr.split(";")) {
                auto items = lineStr.split(",");
                if (items.length != 7) {
                    string errMsg = "The 'extract-line' string requires exactly 7 values.\n";
                    errMsg ~= format("You have provided %d items.\n", items.length);
                    errMsg ~= format("The problematic string is: %s\n", lineStr);
                    throw new Error(errMsg);
                }
                Vector3 p0 = Vector3(to!number(items[0]), to!number(items[1]), to!number(items[2]));
                Vector3 p1 = Vector3(to!number(items[3]), to!number(items[4]), to!number(items[5]));
                size_t n = to!size_t(items[6]);
                auto count = soln.find_enclosing_cells_along_line(p0, p1, n, cells_found);
                writeln("# Info: Found ", count, " cells from point ", p0, " to point ", p1);
            } // end foreach lineStr
            foreach(i; 0 .. cells_found.length) {
                size_t ib = cells_found[i][0]; size_t idx = cells_found[i][1];
                outFile.writeln(soln.flowBlocks[ib].values_as_string(idx));
            }
        } // end foreach tindx
    } // end if extractLineStr

    if (extractSolidLineStr.length > 0) {
        writeln("Extracting data along a straight line between end points in solid domains.");
        // The output may go to a user-specified file, or stdout.
        File outFile;
        if (outputFileName.length > 0) {
            outFile = File(outputFileName, "w");
            writeln("Output will be sent to File: ", outputFileName);
        } else {
            outFile = stdout;
        }
        bool header_written = false;
        extractSolidLineStr = extractSolidLineStr.strip();
        extractSolidLineStr = extractSolidLineStr.replaceAll(regex("\""), "");
        foreach (tindx; tindx_list_to_plot) {
            writeln("  tindx= ", tindx);
            auto soln = new SolidSolution(jobName, ".", tindx, GlobalConfig.nSolidBlocks);
            if (luaRefSoln.length > 0) soln.subtract_ref_soln(luaRefSoln);
            if (!header_written) {
                // Gnuplot column labels.
                outFile.writeln(soln.solidBlocks[0].variable_names_as_string(true));
                header_written = true;
            } else {
                // Gnuplot 2 blank lines separates datasets.
                outFile.write("\n\n");
            }
            size_t[2][] cells_found; // accumulate the identies of the cells found here
            foreach(lineStr; extractSolidLineStr.split(";")) {
                auto items = lineStr.split(",");
                if (items.length != 7) {
                    string errMsg = "The 'extract-solid-line' string requires exactly 7 values.\n";
                    errMsg ~= format("You have provided %d items.\n", items.length);
                    errMsg ~= format("The problematic string is: %s\n", lineStr);
                    throw new Error(errMsg);
                }
                Vector3 p0 = Vector3(to!number(items[0]), to!number(items[1]), to!number(items[2]));
                Vector3 p1 = Vector3(to!number(items[3]), to!number(items[4]), to!number(items[5]));
                size_t n = to!size_t(items[6]);
                auto count = soln.find_enclosing_cells_along_line(p0, p1, n, cells_found);
                writeln("# Info: Found ", count, " cells from point ", p0, " to point ", p1);
            } // end foreach lineStr
            foreach(i; 0 .. cells_found.length) {
                size_t ib = cells_found[i][0]; size_t idx = cells_found[i][1];
                outFile.writeln(soln.solidBlocks[ib].values_as_string(idx));
            }
        } // end foreach tindx
    } // end if extractSolidLineStr
    //
    if (computeLoadsOnGroupStr.length > 0) {
        writeln("Computing loads on group: " ~ computeLoadsOnGroupStr ~ " .");
        // The output may go to a user-specified file, or stdout.
        File outFile;
        if (outputFileName.length > 0) {
            outFile = File(outputFileName, "w");
            writeln("Output will be sent to File: ", outputFileName);
        } else {
            outFile = stdout;
        }
        string groupTag = computeLoadsOnGroupStr;
        foreach (tindx; tindx_list_to_plot) {
            writeln("  tindx= ", tindx);
            auto soln = new FlowSolution(jobName, ".", tindx, GlobalConfig.nFluidBlocks, -1, GlobalConfig.flow_format, tag);
            soln.add_aux_variables(addVarsList, tag);
            if (luaRefSoln.length > 0) soln.subtract_ref_soln(luaRefSoln);
            number Fx = 0.0; number Fy = 0.0; number Fz = 0.0; number F = 0.0; number q = 0.0;
            foreach (blk_indx; 0..GlobalConfig.nFluidBlocks) {
                foreach (boundary_indx; 0..soln.flowBlocks[blk_indx].bcGroups.length) {
                    auto surf_grid = soln.gridBlocks[blk_indx].get_boundary_grid(boundary_indx);
                    size_t[] surf_cells = soln.gridBlocks[blk_indx].get_list_of_boundary_cells(boundary_indx);
                    size_t new_dimensions = surf_grid.dimensions;
                    // The following should work for both structured and unstructured grids.
                    size_t new_nic = max(surf_grid.niv-1, 1);
                    size_t new_njc = max(surf_grid.njv-1, 1);
                    size_t new_nkc = max(surf_grid.nkv-1, 1);
                    assert(new_nic*new_njc*new_nkc == surf_cells.length, "mismatch is number of cells");
                    auto surf_flow = new FluidBlockLite(soln.flowBlocks[blk_indx], surf_cells,
                                                        new_dimensions, new_nic, new_njc, new_nkc);
                    // At this stage we should have a surface flow structure, and a sufrace grid.
                }
            }
        }
    } // end if computeLoadsOnGroupStr
    //
    if (normsStr.length > 0) {
        writeln("Norms for variables.");
        normsStr = normsStr.strip();
        normsStr = normsStr.replaceAll(regex("\""), "");
        foreach (tindx; tindx_list_to_plot) {
            writeln("  tindx= ", tindx);
            auto soln = new FlowSolution(jobName, ".", tindx, GlobalConfig.nFluidBlocks, -1, GlobalConfig.flow_format, tag);
            soln.add_aux_variables(addVarsList, tag);
            if (luaRefSoln.length > 0) soln.subtract_ref_soln(luaRefSoln);
            //
            SolidSolution solidSoln;
            if ( GlobalConfig.nSolidBlocks > 0 ) {
                solidSoln = new SolidSolution(jobName, ".", tindx, GlobalConfig.nSolidBlocks);
                if (luaRefSoln.length > 0) solidSoln.subtract_ref_soln(luaRefSoln);
            }
            //
            // Work on flow blocks first
            writeln("normsStr= ", normsStr);
            foreach (varName; normsStr.split(",")) {
                writeln("flow: varName= ", varName);
                if (!canFind(soln.flowBlocks[0].variableNames, varName)) {
                    writeln(format("Requested variable name \"%s\" not in list of flow variables.", varName));
                    continue;
                }
                auto norms = soln.compute_volume_weighted_norms(varName, regionStr);
                write("    variable= ", varName, "\n");
                write(format(" L1= %.18e L2= %.18e Linf= %.18e\n",
                             norms[0], norms[1], norms[2]));
                write(" x= ", norms[3], " y= ", norms[4], " z= ", norms[5]);
                write("\n");
            } // end foreach varName
            // Then work on solid blocks
            if ( GlobalConfig.nSolidBlocks > 0 ) {
                writeln("normsStr= ", normsStr);
                foreach (varName; normsStr.split(",")) {
                    writeln("solid: varName= ", varName);
                    if (!canFind(solidSoln.solidBlocks[0].variableNames, varName)) {
                        writeln(format("Requested variable name \"%s\" not in list of solid variables.", varName));
                        continue;
                    }
                    auto norms = solidSoln.compute_volume_weighted_norms(varName, regionStr);
                    write("    variable= ", varName, "\n");
                    write(format(" L1= %.18e L2= %.18e Linf= %.18e\n",
                                 norms[0], norms[1], norms[2]));
                    write(" x= ", norms[3], " y= ", norms[4], " z= ", norms[5]);
                    write("\n");
                } // end foreach varName
            } // end if nSolidBlocks > 0
        } // end foreach tindx
    } // end if normsStr
    //
} // end post_process()

//-----------------------------------------------------------------------

size_t[] decode_range_indices(string rangeStr, size_t first, size_t endplus1)
// Decode strings such as "0:$", ":", "0:3", "$"
// On input, first and endplus1 represent the largest, available range.
// Return the pair of numbers that can be used in a foreach loop range.
{
    if (rangeStr == ":") {
        return [first, endplus1];
    }
    if (canFind(rangeStr, ":")) {
        // We have a range specification to pull apart.
        auto items = rangeStr.split(":");
        first = to!size_t(items[0]);
        if (items.length > 1 && items[1] != "$") {
            // Presume that we have a second integer.
            size_t new_endplus1 = to!size_t(items[1]);
            if (new_endplus1 < endplus1) endplus1 = new_endplus1;
        }
    } else if (rangeStr == "$") {
        // Wit just a single "$" specified, we want only the last index.
        first = endplus1 - 1;
    }else {
        // Presume that we have a single integer.
        first = to!size_t(rangeStr);
        if (first < endplus1) endplus1 = first+1;
    }
    return [first, endplus1];
} // end decode_range_indices()

double[int] readTimesFile(string jobName, string jobDir=".")
{
    double[int] times_dict;
    // Read the times file for all tindx values.
    auto timesFile = File(jobDir ~ "/config/" ~ jobName ~ ".times");
    auto line = timesFile.readln().strip();
    while (line.length > 0) {
        if (line[0] != '#') {
            // Process a non-comment line.
            auto tokens = line.split();
            times_dict[to!int(tokens[0])] = to!double(tokens[1]);
        }
        line = timesFile.readln().strip();
    }
    timesFile.close();
    return times_dict;
} // end readTimesFile()
