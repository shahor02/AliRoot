// LesHouches.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the LHAup and
// LHAupLHEF classes.

#include "Pythia8/LesHouches.h"

// Access time information.
#include <ctime>

namespace Pythia8 {

//==========================================================================

// LHAup class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// LHA convention with cross section in pb may require conversion from mb.
const double LHAup::CONVERTMB2PB = 1e9;

//--------------------------------------------------------------------------

// Print the initialization info; to check it worked.

void LHAup::listInit(ostream& os) {

  // Header.
  os << "\n --------  LHA initialization information  ------------ \n";

  // Beam info.
  os << fixed << setprecision(3)
     << "\n  beam    kind      energy  pdfgrp  pdfset \n"
     << "     A  " << setw(6) << idBeamASave
     <<  setw(12) << eBeamASave
     << setw(8) << pdfGroupBeamASave
     << setw(8) << pdfSetBeamASave << "\n"
     << "     B  " << setw(6) << idBeamBSave
     <<  setw(12) << eBeamBSave
     << setw(8) << pdfGroupBeamBSave
     << setw(8) << pdfSetBeamBSave << "\n";

  // Event weighting strategy.
  os << "\n  Event weighting strategy = " << setw(2)
     << strategySave << "\n" ;

  // Process list.
  os << scientific << setprecision(4)
     << "\n  Processes, with strategy-dependent cross section info \n"
     << "  number      xsec (pb)      xerr (pb)      xmax (pb) \n" ;
  for (int ip = 0; ip < int(processes.size()); ++ip) {
    os << setw(8) << processes[ip].idProc
       << setw(15) << processes[ip].xSecProc
       << setw(15) << processes[ip].xErrProc
       << setw(15) << processes[ip].xMaxProc << "\n";
  }

  // Finished.
  os << "\n --------  End LHA initialization information  -------- \n";

}

//--------------------------------------------------------------------------

// Print the event info; to check it worked.

void LHAup::listEvent(ostream& os) {

  // Header.
  os << "\n --------  LHA event information and listing  -------------"
     << "--------------------------------------------------------- \n";

  // Basic event info.
  os << scientific << setprecision(4)
     << "\n    process = " << setw(8) << idProc
     << "    weight = " << setw(12) << weightProc
     << "     scale = " << setw(12) << scaleProc << " (GeV) \n"
     << "                   "
     << "     alpha_em = " << setw(12) << alphaQEDProc
     << "    alpha_strong = " << setw(12) << alphaQCDProc << "\n";

  // Particle list
  os << fixed << setprecision(3)
     << "\n    Participating Particles \n"
     << "    no        id stat     mothers     colours      p_x        "
     << "p_y        p_z         e          m        tau    spin \n" ;
  for (int ip = 1; ip < int(particles.size()); ++ip) {
    os << setw(6) << ip
       << setw(10) << particles[ip].idPart
       << setw(5) << particles[ip].statusPart
       << setw(6) << particles[ip].mother1Part
       << setw(6) << particles[ip].mother2Part
       << setw(6) << particles[ip].col1Part
       << setw(6) << particles[ip].col2Part
       << setw(11) << particles[ip].pxPart
       << setw(11) << particles[ip].pyPart
       << setw(11) << particles[ip].pzPart
       << setw(11) << particles[ip].ePart
       << setw(11) <<  particles[ip].mPart
       << setw(8) <<  particles[ip].tauPart
       << setw(8) <<  particles[ip].spinPart << "\n";
  }

  // PDF info - optional.
  if (pdfIsSetSave) os << "\n     pdf: id1 =" << setw(5) << id1pdfSave
    << " id2 =" << setw(5) << id2pdfSave
    << " x1 ="  << scientific << setw(10) << x1pdfSave
    << " x2 =" << setw(10) << x2pdfSave
    << " scalePDF =" << setw(10) << scalePDFSave
    << " pdf1 =" << setw(10) << pdf1Save
    << " pdf2 =" << setw(10) << pdf2Save << "\n";

  // Finished.
  os << "\n --------  End LHA event information and listing  ---------"
     << "--------------------------------------------------------- \n";

}

//--------------------------------------------------------------------------

// Open and write header to a Les Houches Event File.

bool LHAup::openLHEF(string fileNameIn) {

  // Open file for writing. Reset it to be empty.
  fileName = fileNameIn;
  const char* cstring = fileName.c_str();
  osLHEF.open(cstring, ios::out | ios::trunc);
  if (!osLHEF) {
    infoPtr->errorMsg("Error in LHAup::openLHEF:"
      " could not open file", fileName);
    return false;
  }

  // Read out current date and time.
  time_t t = time(0);
  strftime(dateNow,12,"%d %b %Y",localtime(&t));
  strftime(timeNow,9,"%H:%M:%S",localtime(&t));

  // Write header.
  osLHEF << "<LesHouchesEvents version=\"1.0\">\n"
         << "<!--\n"
         << "  File written by Pythia8::LHAup on "
         << dateNow << " at " << timeNow << "\n"
         << "-->" << endl;

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Write initialization information to a Les Houches Event File.

bool LHAup::initLHEF() {

  // Write information on beams.
  osLHEF << "<init>\n" << scientific << setprecision(6)
         << "  " << idBeamASave       << "  " << idBeamBSave
         << "  " << eBeamASave        << "  " << eBeamBSave
         << "  " << pdfGroupBeamASave << "  " << pdfGroupBeamBSave
         << "  " << pdfSetBeamASave   << "  " << pdfSetBeamBSave
         << "  " << strategySave      << "  " << processes.size() << "\n";

  // Write information on all the subprocesses.
  for (int ip = 0; ip < int(processes.size()); ++ip)
    osLHEF << " " << setw(13) << processes[ip].xSecProc
           << " " << setw(13) << processes[ip].xErrProc
           << " " << setw(13) << processes[ip].xMaxProc
           << " " << setw(6) << processes[ip].idProc << "\n";

  // Done.
  osLHEF << "</init>" << endl;
  return true;

}

//--------------------------------------------------------------------------

// Write event information to a Les Houches Event File.
// Normal mode is to line up event info in columns, but the non-verbose
// altnernative saves space at the expense of human readability.

bool LHAup::eventLHEF(bool verbose) {

  // Default verbose option.
  if (verbose) {

    // Write information on process as such.
    osLHEF << "<event>\n" << scientific << setprecision(6)
           << " " << setw(5) << particles.size() - 1
           << " " << setw(5) << idProc
           << " " << setw(13) << weightProc
           << " " << setw(13) << scaleProc
           << " " << setw(13) << alphaQEDProc
           << " " << setw(13) << alphaQCDProc << "\n";

    // Write information on the particles, excluding zeroth.
    for (int ip = 1; ip < int(particles.size()); ++ip) {
      LHAParticle& ptNow = particles[ip];
      osLHEF << " " << setw(8) << ptNow.idPart
             << " " << setw(5) << ptNow.statusPart
             << " " << setw(5) << ptNow.mother1Part
             << " " << setw(5) << ptNow.mother2Part
             << " " << setw(5) << ptNow.col1Part
             << " " << setw(5) << ptNow.col2Part << setprecision(10)
             << " " << setw(17) << ptNow.pxPart
             << " " << setw(17) << ptNow.pyPart
             << " " << setw(17) << ptNow.pzPart
             << " " << setw(17) << ptNow.ePart
             << " " << setw(17) <<  ptNow.mPart << setprecision(6);
      if (ptNow.tauPart == 0.) osLHEF << " 0.";
      else osLHEF << " " << setw(13) << ptNow.tauPart;
      if (ptNow.spinPart == 9.) osLHEF << " 9.";
      else osLHEF << " " << setw(13) << ptNow.spinPart;
      osLHEF << "\n";
    }

    // Optionally write information on PDF values at hard interaction.
    if (pdfIsSetSave) osLHEF << "#pdf"
             << " " << setw(4) << id1pdfSave
             << " " << setw(4) << id2pdfSave
             << " " << setw(13) << x1pdfSave
             << " " << setw(13) << x2pdfSave
             << " " << setw(13) << scalePDFSave
             << " " << setw(13) << pdf1Save
             << " " << setw(13) << pdf2Save << "\n";

  // Alternative non-verbose option.
  } else {

    // Write information on process as such.
    osLHEF << "<event>\n" << scientific << setprecision(6)
           << particles.size() - 1 << " " << idProc       << " "
           << weightProc           << " " << scaleProc    << " "
           << alphaQEDProc         << " " << alphaQCDProc << "\n";

    // Write information on the particles, excluding zeroth.
    for (int ip = 1; ip < int(particles.size()); ++ip) {
      LHAParticle& ptNow = particles[ip];
      osLHEF        << ptNow.idPart      << " " << ptNow.statusPart
             << " " << ptNow.mother1Part << " " << ptNow.mother2Part
             << " " << ptNow.col1Part    << " " << ptNow.col2Part
             << setprecision(10)         << " " << ptNow.pxPart
             << " " << ptNow.pyPart      << " " << ptNow.pzPart
             << " " << ptNow.ePart       << " " << ptNow.mPart
             << setprecision(6);
      if (ptNow.tauPart == 0.) osLHEF << " 0.";
      else osLHEF << " " << setw(13) << ptNow.tauPart;
      if (ptNow.spinPart == 9.) osLHEF << " 9.";
      else osLHEF << " " << setw(13) << ptNow.spinPart;
      osLHEF << "\n";
    }

    // Optionally write information on PDF values at hard interaction.
    if (pdfIsSetSave) osLHEF << "#pdf" << " " << id1pdfSave
             << " " << id2pdfSave << " " << x1pdfSave << " " << x2pdfSave
             << " " << scalePDFSave << " " << pdf1Save << " " << pdf2Save
             << "\n";
  }

  // Done.
  osLHEF << "</event>" << endl;
  return true;

}

//--------------------------------------------------------------------------

// Write end of a Les Houches Event File and close it.

bool LHAup::closeLHEF(bool updateInit) {

  // Write an end to the file.
  osLHEF << "</LesHouchesEvents>" << endl;
  osLHEF.close();

  // Optionally update the cross section information.
  if (updateInit) {
    const char* cstring = fileName.c_str();
    osLHEF.open(cstring, ios::in | ios::out);

    // Rewrite header; identically with what openLHEF did.
    osLHEF << "<LesHouchesEvents version=\"1.0\">\n"
           << "<!--\n"
           << "  File written by Pythia8::LHAup on "
           << dateNow << " at " << timeNow << "\n"
           << "-->" << endl;

    // Redo initialization information.
    initLHEF();
    osLHEF.close();
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Read in initialization information from a Les Houches Event File.

bool LHAup::setInitLHEF(istream& is, bool readHeaders) {

  // Check that first line is consistent with proper LHEF file.
  string line;
  if (!getline(is, line)) return false;
  if (line.find("<LesHouchesEvents") == string::npos) return false;
  if (line.find("version=\"1.0\"" ) == string::npos ) return false;

  // What to search for if reading headers; if not reading
  // headers then return to default behaviour
  string headerTag = (readHeaders) ? "<header>" : "<init";

  // Loop over lines until an <init (or optionally <header>) tag
  // is found first on a line.
  string tag = " ";
  do {
    if (!getline(is, line)) return false;
    if (line.find_first_not_of(" \n\t\v\b\r\f\a") != string::npos) {
      istringstream getfirst(line);
      getfirst >> tag;
      if (!getfirst) return false;
    }
  } while (tag != "<init>" && tag != "<init" && tag != headerTag);

  // If header tag found, process if required
  if (readHeaders == true && tag == headerTag) {
    // Temporary local storage
    map < string, string > headerMap;

    // Loop over lines until an <init> tag is found.
    bool read = true, newKey = false;
    string key = "base";
    vector < string > keyVec;

    while (true) {
      if (!getline(is, line)) return false;

      // Break lines containing multiple tags into two segments.
      // (Could be generalized to multiple segments but this is
      // sufficient to handle at least <tag>info</tag> on same line.
      size_t firstTagEnd    = line.find_first_of(">");
      size_t secondTagBegin = line.find_first_of("<",firstTagEnd);
      vector<string> lineVec;
      if (firstTagEnd != string::npos && secondTagBegin != string::npos) {
        lineVec.push_back(line.substr(0,secondTagBegin));
        lineVec.push_back(line.substr(secondTagBegin,
          line.size()-secondTagBegin));
      }
      else {
        lineVec.push_back(line);
      }

      // Loop over segments of current line
      for (int iVec=0; iVec<int(lineVec.size()); ++iVec) {
        line = lineVec[iVec];

        // Clean line to contain only valid characters
        size_t posBeg = line.find_first_not_of(" \n\t\v\b\r\f\a");
        size_t posEnd = line.find_last_not_of(" \n\t\v\b\r\f\a");
        string lineClean = " ";
        if (posBeg != string::npos && posEnd != string::npos && posBeg
          < posEnd) {
          lineClean = line.substr(posBeg, posEnd - posBeg + 1);
          posBeg = 0;
          posEnd = lineClean.size();
        }

        // Check for empty line
        if (lineClean == " " || posBeg >= posEnd) continue;

        // PZS Jan 2015: Allow multiple open/close tags on a single line.
        size_t tagBeg =  lineClean.find_first_of("<");
        size_t tagEnd =  lineClean.find_first_of(">");

        while (tagBeg != string::npos && tagBeg < tagEnd) {

          // Update remainder (non-tag) part of line, for later storage
          posBeg = tagEnd+1;

          // Only take the first word of the tag,
          tag = lineClean.substr(tagBeg + 1, tagEnd - tagBeg - 1);
          istringstream getfirst(tag);
          getfirst >> tag;

          // Prepare for next while iteration:
          // Look for next tag on line and update posBeg and posEnd.
          tagBeg = lineClean.find_first_of("<",tagEnd);
          tagEnd = lineClean.find_first_of(">",tagBeg+1);

          // Tag present, so handle here
          if (getfirst) {

            // Exit condition
            if (tag == "init") break;

            // End of header block; keep reading until <init> tag,
            // but do not store any further information
            else if (tag == "/header") {
              read = false;
              continue;

              // Opening tag
            } else if (tag[0] != '/') {
              keyVec.push_back(tag);
              newKey = true;
              continue;

              // Closing tag that matches current key
            } else if (tag == "/" + keyVec.back()) {
              keyVec.pop_back();
              newKey = true;
              continue;

              // Also check for forgotten close tag: next-to-last element
            } else if (keyVec.size() >= 2
                       && tag == "/" + keyVec[keyVec.size()-2]) {
              infoPtr->errorMsg("Warning in LHAup::setInitLHEF:"
                                " corrupt LHEF end tag",keyVec.back());
              keyVec.pop_back();
              keyVec.pop_back();
              newKey = true;
              continue;
            }

          } // if (getfirst)

        } // Loop over tags

        // Exit condition
        if (tag == "init") break;

        // At this point the (rest of) the line is not a tag;
        // If no longer reading anything, skip.
        if (!read) continue;

        // Check for key change
        if (newKey) {
          if (keyVec.empty()) key = "base";
          else                key = keyVec[0];
          for (size_t i = 1; i < keyVec.size(); i++)
            key += "." + keyVec[i];
          newKey = false;
        }

        // Check if anything remains to store of this line
        posBeg = line.find_first_not_of(" \n\t\v\b\r\f\a",posBeg);
        if (posBeg == string::npos || posBeg > posEnd) continue;

        // Append information to local storage
        headerMap[key] += line.substr(posBeg,posEnd - posBeg + 1) + "\n";

      } // Loop over line segments

      // Exit condition
      if (tag == "init") break;

    } // while (true)

    // Copy information to info using LHAup::setInfoHeader
    for (map < string, string >::iterator it = headerMap.begin();
         it != headerMap.end(); it++)
      setInfoHeader(it->first, it->second);

  } // if (readHeaders == true && tag == headerTag)

  // Read in first info line; done if empty.
  if (!getline(is, line)) return false;
  if (line.find("</init") != string::npos) return true;

  // Read in beam and strategy info, and store it.
  int idbmupA, idbmupB;
  double ebmupA, ebmupB;
  int pdfgupA, pdfgupB, pdfsupA, pdfsupB, idwtup, nprup;
  istringstream getbms(line);
  getbms >> idbmupA >> idbmupB >> ebmupA >> ebmupB >> pdfgupA
     >> pdfgupB >> pdfsupA >> pdfsupB >> idwtup >> nprup;
  if (!getbms) return false;
  setBeamA(idbmupA, ebmupA, pdfgupA, pdfsupA);
  setBeamB(idbmupB, ebmupB, pdfgupB, pdfsupB);
  setStrategy(idwtup);

  // Read in process info, one process at a time, and store it.
  double xsecup, xerrup, xmaxup;
  xSecSumSave = 0.;
  xErrSumSave = 0.;
  int lprup;
  for (int ip = 0; ip < nprup; ++ip) {
    if (!getline(is, line)) return false;
    istringstream getpro(line);
    getpro >> xsecup >> xerrup >> xmaxup >> lprup ;
    if (!getpro) return false;
    addProcess(lprup, xsecup, xerrup, xmaxup);
    xSecSumSave += xsecup;
    xErrSumSave += pow2(xerrup);
  }
  xErrSumSave = sqrt(xErrSumSave);

  // Reading worked.
  return true;

}

//--------------------------------------------------------------------------

// Read in event information from a Les Houches Event File,
// into a staging area where it can be reused by setOldEventLHEF.

bool LHAup::setNewEventLHEF(istream& is) {

  // Loop over lines until an <event tag is found first on a line.
  string line, tag;
  do {
    if (!getline(is, line)) return false;
    if (line.find_first_not_of(" \n\t\v\b\r\f\a") != string::npos) {
      istringstream getfirst(line);
      getfirst >> tag;
      if (!getfirst) return false;
    }
  } while (tag != "<event>" && tag != "<event");

  // Read in process info and store it.
  if (!getline(is, line)) return false;
  istringstream getpro(line);
  getpro >> nupSave >> idprupSave >> xwgtupSave >> scalupSave
    >> aqedupSave >> aqcdupSave;
  if (!getpro) return false;

  // Reset particlesSave vector, add slot-0 empty particle.
  particlesSave.clear();
  particlesSave.push_back( LHAParticle() );

  // Read in particle info one by one, and store it.
  // Note unusual C++ loop range, to better reflect LHA/Fortran standard.
  // (Recall that process(...) above added empty particle at index 0.)
  int idup, istup, mothup1, mothup2, icolup1, icolup2;
  double pup1, pup2, pup3, pup4, pup5, vtimup, spinup;
  for (int ip = 1; ip <= nupSave; ++ip) {
    if (!getline(is, line)) return false;
    istringstream getall(line);
    getall >> idup >> istup >> mothup1 >> mothup2 >> icolup1 >> icolup2
      >> pup1 >> pup2 >> pup3 >> pup4 >> pup5 >> vtimup >> spinup;
    if (!getall) return false;
    particlesSave.push_back( LHAParticle( idup, istup, mothup1, mothup2,
      icolup1, icolup2, pup1, pup2, pup3, pup4, pup5, vtimup, spinup, -1.) );
  }

  // Flavour and x values of hard-process initiators.
  id1InSave = particlesSave[1].idPart;
  id2InSave = particlesSave[2].idPart;
  x1InSave  = (eBeamASave > 0.) ? particlesSave[1].ePart / eBeamASave : 0.;
  x2InSave  = (eBeamBSave > 0.) ? particlesSave[2].ePart / eBeamBSave : 0.;

  // Continue parsing till </event>. Look for optional info on the way.
  getPDFSave = false;
  getScale   = false;
  do {
    if (!getline(is, line)) return false;
    istringstream getinfo(line);
    getinfo >> tag;
    if (!getinfo) return false;

    // Extract PDF info if present.
    if (tag == "#pdf" && !getPDFSave) {
      getinfo >> id1pdfInSave >> id2pdfInSave >> x1pdfInSave >> x2pdfInSave
              >> scalePDFInSave >> pdf1InSave >> pdf2InSave;
      if (!getinfo) return false;
      getPDFSave = true;

    // Extract scale info if present.
    } else if (tag == "#" && !getScale) {
      double scaleIn = 0;
      for (int i = 3; i < int(particlesSave.size()); ++i)
      if (particlesSave[i].statusPart == 1) {
        if ( !(getinfo >> scaleIn) ) return false;
        particlesSave[i].scalePart = scaleIn;
      }
      if (!getinfo) return false;
      getScale = true;
    }
  } while (tag != "</event>" && tag != "</event");

  // Need id and x values even when no PDF info. Rest empty.
  if (!getPDFSave) {
    id1pdfInSave   = id1InSave;
    id2pdfInSave   = id2InSave;
    x1pdfInSave    = x1InSave;
    x2pdfInSave    = x2InSave;
    scalePDFInSave = 0.;
    pdf1InSave     = 0.;
    pdf2InSave     = 0.;
  }

  // Reading worked.
  return true;

}

//--------------------------------------------------------------------------

// Make current event information read in by setNewEventLHEF.

bool LHAup::setOldEventLHEF() {

  // Store saved event, optionally also parton density information.
  setProcess( idprupSave, xwgtupSave, scalupSave, aqedupSave, aqcdupSave);
  for (int ip = 1; ip <= nupSave; ++ip) addParticle( particlesSave[ip] );
  setIdX( id1InSave, id2InSave, x1InSave, x2InSave);
  setPdf( id1pdfInSave, id2pdfInSave, x1pdfInSave, x2pdfInSave,
    scalePDFInSave, pdf1InSave, pdf2InSave, getPDFSave);

  // Done;
  return true;

}

//--------------------------------------------------------------------------

// Open a file using provided ifstream and return a pointer to an istream
// that can be used to process the file.

istream* LHAup::openFile(const char *fn, ifstream &ifs) {
  // Open the file
  ifs.open(fn);
  return (istream *) &ifs;
}

//--------------------------------------------------------------------------

// Companion method to 'openFile', above.
// Correctly deallocates memory if required before closing the file.

void LHAup::closeFile(istream *&is, ifstream &ifs) {
  // If the istream pointer is not NULL and is not the
  // same as the ifstream, then delete pointer.
  if (is && is != &ifs) delete is;
  is = NULL;

  // Close the file
  if (ifs.is_open()) ifs.close();
}

//==========================================================================

// LHAupLHEF class.

//--------------------------------------------------------------------------

// Routine for doing the job of reading and setting initialization info.

bool LHAupLHEF::setInitLHEF( istream & isIn, bool readHead ) {

  // Done if there was a problem with initialising the reader
  if (!reader.isGood) return false;

  // Construct header information (stored in comments strings or optional
  // header file), so that reading of headers is possible.
  string comments;
  comments+="<LesHouchesEvents version =\"3.0\">\n";
  comments+="<header>\n";
  comments+=reader.headerComments;
  comments+="</header>\n";
  comments+="<init>\n";
  comments+=reader.initComments;
  comments+="</init>\n";
  istringstream is1(comments);
  bool useComments = (headerfile == NULL);
  istream & iss((useComments ? is1 : isIn));

  // Check that first line is consistent with proper LHEF file.
  string line;
  if ( useComments && !getline(iss,line)) return false;
  if (!useComments && !getLine(line)) return false;

  // What to search for if reading headers; if not reading
  // headers then return to default behaviour
  string headerTag = (readHead) ? "<header>" : "<init";

  // Loop over lines until an <init (or optionally <header>) tag
  // is found first on a line.
  string tag = " ";
  do {
    if ( useComments && !getline(iss,line)) return false;
    if (!useComments && !getLine(line)) return false;
    if (line.find_first_not_of(" \n\t\v\b\r\f\a") != string::npos) {
      istringstream getfirst(line);
      getfirst >> tag;
      if (!getfirst) return false;
    }
  } while (tag != "<init>" && tag != "<init" && tag != headerTag);

  // If header tag found, process if required
  if (readHead == true && tag == headerTag) {
    // Temporary local storage
    map < string, string > headerMap;

    // Loop over lines until an <init> tag is found.
    bool read = true, newKey = false;
    int commentDepth = 0;
    string key = "base";
    vector < string > keyVec;
    while (true) {
      if ( useComments && !getline(iss,line)) return false;
      if (!useComments && !getLine(line)) return false;

      // Tell XML parser to ignore comment and CDATA blocks
      // If we are currently inside a comment block, check for block end
      if (commentDepth >= 1 && line.find("-->") != string::npos)
        commentDepth--;
      if (commentDepth >= 1 && line.find("]]>") != string::npos)
        commentDepth--;
      // If the comment block did not end on this line, skip to next line
      if (commentDepth >= 1) continue;
      // Check for beginning of comment blocks (parse until comment begins)
      if (line.find("<!--") != string::npos) {
        if (line.find("-->") == string::npos) commentDepth++;
        int comBeg = line.find("<!--");
        line = line.substr(0,comBeg);
      }
      // Check for beginning of CDATA statement  (parse until CDATA begins)
      if (line.find("<![cdata[") != string::npos
          || line.find("<![CDATA[") != string::npos) {
        if (line.find("]]>") == string::npos) commentDepth++;
        int comBeg = line.find("<![");
        line = line.substr(0,comBeg);
      }

      // Break lines containing multiple tags into two segments.
      // (Could be generalized to multiple segments but this is
      // sufficient to handle at least <tag>info</tag> on same line.
      size_t firstTagEnd    = line.find_first_of(">");
      size_t secondTagBegin = line.find_first_of("<",firstTagEnd);
      vector<string> lineVec;
      if (firstTagEnd != string::npos && secondTagBegin != string::npos) {
        lineVec.push_back(line.substr(0,secondTagBegin));
        lineVec.push_back(line.substr(secondTagBegin,
          line.size()-secondTagBegin));
      }
      else {
        lineVec.push_back(line);
      }

      // Loop over segments of current line
      for (int iVec=0; iVec<int(lineVec.size()); ++iVec) {
        line = lineVec[iVec];

        // Clean line to contain only valid characters
        size_t posBeg = line.find_first_not_of(" \n\t\v\b\r\f\a");
        size_t posEnd = line.find_last_not_of(" \n\t\v\b\r\f\a");
        string lineClean = " ";
        if (posBeg != string::npos && posEnd != string::npos && posBeg
          < posEnd) {
          lineClean = line.substr(posBeg, posEnd - posBeg + 1);
          posBeg = 0;
          posEnd = lineClean.size();
        }

        // Check for empty line
        if (lineClean == " " || posBeg >= posEnd) continue;

        // PZS Jan 2015: Allow multiple open/close tags on a single line.
        size_t tagBeg =  lineClean.find_first_of("<");
        size_t tagEnd =  lineClean.find_first_of(">");

        while (tagBeg != string::npos && tagBeg < tagEnd) {

          // Update remainder (non-tag) part of line, for later storage
          posBeg = tagEnd+1;

          // Only take the first word of the tag,
          tag = lineClean.substr(tagBeg + 1, tagEnd - tagBeg - 1);
          istringstream getfirst(tag);
          getfirst >> tag;

          // Prepare for next while iteration:
          // Look for next tag on line and update posBeg and posEnd.
          tagBeg = lineClean.find_first_of("<",tagEnd);
          tagEnd = lineClean.find_first_of(">",tagBeg+1);

          // Tag present, so handle here
          if (getfirst) {

            // Exit condition
            if (tag == "init") break;

            // End of header block; keep reading until <init> tag,
            // but do not store any further information
            else if (tag == "/header") {
              read = false;
              continue;

              // Opening tag
            } else if (tag[0] != '/') {
              keyVec.push_back(tag);
              newKey = true;
              continue;

              // Closing tag that matches current key
            } else if (tag == "/" + keyVec.back()) {
              keyVec.pop_back();
              newKey = true;
              continue;

              // Also check for forgotten close tag: next-to-last element
            } else if (keyVec.size() >= 2
                       && tag == "/" + keyVec[keyVec.size()-2]) {
              infoPtr->errorMsg("Warning in LHAupLHEF::setInitLHEF:"
                                " corrupt LHEF end tag",keyVec.back());
              keyVec.pop_back();
              keyVec.pop_back();
              newKey = true;
              continue;
            }

          } // if (getfirst)

        } // Loop over tags

        // Exit condition
        if (tag == "init") break;

        // At this point the (rest of) the line is not a tag;
        // If no longer reading anything, skip.
        if (!read) continue;

        // Check for key change
        if (newKey) {
          if (keyVec.empty()) key = "base";
          else                key = keyVec[0];
          for (size_t i = 1; i < keyVec.size(); i++)
            key += "." + keyVec[i];
          newKey = false;
        }

        // Check if anything remains to store of this line
        posBeg = line.find_first_not_of(" \n\t\v\b\r\f\a",posBeg);
        if (posBeg == string::npos || posBeg > posEnd) continue;

        // Append information to local storage
        headerMap[key] += line.substr(posBeg,posEnd - posBeg + 1) + "\n";

      } // Loop over line segments

      // Exit condition
      if (tag == "init") break;

    } // while (true)

    // Copy information to info using LHAup::setInfoHeader
    for (map < string, string >::iterator it = headerMap.begin();
         it != headerMap.end(); it++)
      setInfoHeader(it->first, it->second);

  } // if (readHead == true && tag == headerTag)

  // Extract beam and strategy info, and store it.
  int idbmupA, idbmupB;
  double ebmupA, ebmupB;
  int pdfgupA, pdfgupB, pdfsupA, pdfsupB, idwtup, nprup;

  idbmupA  = reader.heprup.IDBMUP.first;
  idbmupB  = reader.heprup.IDBMUP.second;
  ebmupA   = reader.heprup.EBMUP.first;
  ebmupB   = reader.heprup.EBMUP.second;
  pdfgupA  = reader.heprup.PDFGUP.first;
  pdfgupB  = reader.heprup.PDFGUP.first;
  pdfsupA  = reader.heprup.PDFSUP.first;
  pdfsupB  = reader.heprup.PDFSUP.second;
  idwtup   = reader.heprup.IDWTUP;
  nprup    = reader.heprup.NPRUP;

  setBeamA(idbmupA, ebmupA, pdfgupA, pdfsupA);
  setBeamB(idbmupB, ebmupB, pdfgupB, pdfsupB);
  setStrategy(idwtup);

  // Read in process info, one process at a time, and store it.
  double xsecup, xerrup, xmaxup;
  xSecSumSave = 0.;
  xErrSumSave = 0.;
  int lprup;
  infoPtr->sigmaLHEFSave.resize(0);
  for (int ip = 0; ip < nprup; ++ip) {
    xsecup = reader.heprup.XSECUP[ip];
    xerrup = reader.heprup.XERRUP[ip];
    xmaxup = reader.heprup.XMAXUP[ip];
    lprup = reader.heprup.LPRUP[ip];
    addProcess(lprup, xsecup, xerrup, xmaxup);
    xSecSumSave += xsecup;
    xErrSumSave += pow(xerrup,2);
    infoPtr->sigmaLHEFSave.push_back(xsecup);
  }
  xErrSumSave = sqrt(xErrSumSave);

  // Now extract LHEF 2.0/3.0 novelties.
  infoPtr->setLHEF3InitInfo();
  if (reader.version > 1) {
    infoPtr->setLHEF3InitInfo( reader.version,
      &reader.heprup.initrwgt, &(reader.heprup.generators),
      &(reader.heprup.weightgroups), &(reader.heprup.weights));
  }

  // Reading worked.
  return true;
}

//--------------------------------------------------------------------------

// Routine for doing the job of reading and setting info on next event.

bool LHAupLHEF::setNewEventLHEF() {

  // Done if the reader finished preemptively.
  if(!reader.readEvent()) return false;

  // Extract process info and store it.
  nupSave     = reader.hepeup.NUP;
  idprupSave  = reader.hepeup.IDPRUP;
  xwgtupSave  = reader.hepeup.XWGTUP;
  scalupSave  = reader.hepeup.SCALUP;
  aqedupSave  = reader.hepeup.AQEDUP;
  aqcdupSave  = reader.hepeup.AQCDUP;

  // Reset particlesSave vector, add slot-0 empty particle.
  particlesSave.clear();
  particlesSave.push_back( Pythia8::LHAParticle() );

  // Read in particle info one by one, and store it.
  // Note unusual C++ loop range, to better reflect LHA/Fortran standard.
  // (Recall that process(...) above added empty particle at index 0.)
  int idup, istup, mothup1, mothup2, icolup1, icolup2;
  double pup1, pup2, pup3, pup4, pup5, vtimup, spinup;
  for ( int i = 0; i < reader.hepeup.NUP; ++i ) {
    // Extract information stored in reader.
    idup     = reader.hepeup.IDUP[i];
    istup    = reader.hepeup.ISTUP[i];
    mothup1  = reader.hepeup.MOTHUP[i].first;
    mothup2  = reader.hepeup.MOTHUP[i].second;
    icolup1  = reader.hepeup.ICOLUP[i].first;
    icolup2  = reader.hepeup.ICOLUP[i].second;
    pup1     = reader.hepeup.PUP[i][0];
    pup2     = reader.hepeup.PUP[i][1];
    pup3     = reader.hepeup.PUP[i][2];
    pup4     = reader.hepeup.PUP[i][3];
    pup5     = reader.hepeup.PUP[i][4];
    vtimup   = reader.hepeup.VTIMUP[i];
    spinup   = reader.hepeup.SPINUP[i];
    particlesSave.push_back( Pythia8::LHAParticle( idup,istup,mothup1,mothup2,
      icolup1, icolup2, pup1, pup2, pup3, pup4, pup5, vtimup, spinup, -1.) );
  }

  // Flavour and x values of hard-process initiators.
  id1InSave = particlesSave[1].idPart;
  id2InSave = particlesSave[2].idPart;
  x1InSave  = (eBeamA() > 0.) ? particlesSave[1].ePart / eBeamA() : 0.;
  x2InSave  = (eBeamB() > 0.) ? particlesSave[2].ePart / eBeamB() : 0.;

  // Parse event comments and look for optional info on the way.
  std::string line, tag;
  std::stringstream ss(reader.eventComments);
  getPDFSave = false;
  getScale   = false;
  getScale   = (setScalesFromLHEF && reader.version == 1) ? false : true;
  while (getline(ss, line)) {
    istringstream getinfo(line);
    getinfo >> tag;
    if (!getinfo) break;
    // Extract PDF info if present.
    if (tag == "#pdf" && !getPDFSave) {
      getinfo >> id1pdfInSave >> id2pdfInSave >> x1pdfInSave >> x2pdfInSave
              >> scalePDFInSave >> pdf1InSave >> pdf2InSave;
      if (!getinfo) return false;
      getPDFSave = true;
    // Extract scale info if present.
    } else if (tag == "#" && !getScale) {
      double scaleIn = 0;
      for (int i = 3; i < int(particlesSave.size()); ++i)
        if (particlesSave[i].statusPart == 1) {
          if ( !(getinfo >> scaleIn) ) return false;
          particlesSave[i].scalePart = scaleIn;
        }
      if (!getinfo) return false;
      getScale = true;
    }
  }

  // Set production scales from <scales> tag.
  if ( setScalesFromLHEF && reader.version > 1 ){
    for ( map<string,double>::const_iterator
      it  = reader.hepeup.scales.attributes.begin();
      it != reader.hepeup.scales.attributes.end(); ++it ) {
      if ( it->first.find_last_of("_") != string::npos) {
        unsigned iFound = it->first.find_last_of("_") + 1;
        int iPos = atoi(it->first.substr(iFound).c_str());
        // Only set production scales of final particles.
        if ( iPos < int(particlesSave.size())
          && particlesSave[iPos].statusPart == 1)
          particlesSave[iPos].scalePart = it->second;
      }
    }
  }

  // Need id and x values even when no PDF info. Rest empty.
  if (!getPDFSave) {
    id1pdfInSave   = id1InSave;
    id2pdfInSave   = id2InSave;
    x1pdfInSave    = x1InSave;
    x2pdfInSave    = x2InSave;
    scalePDFInSave = 0.;
    pdf1InSave     = 0.;
    pdf2InSave     = 0.;
  }

  // Now extract LHEF 2.0/3.0 novelties.
  infoPtr->setLHEF3EventInfo();
  // Set everything for 2.0 and 3.0
  if (reader.version > 1) {
    infoPtr->setLHEF3EventInfo( &reader.hepeup.attributes,
      &reader.hepeup.weights_detailed, &reader.hepeup.weights_compressed,
      &reader.hepeup.scales, &reader.hepeup.weights, &reader.hepeup.rwgt);
  // Try to at least set the event attributes for 1.0
  } else {
    infoPtr->setLHEF3EventInfo( &reader.hepeup.attributes, 0, 0, 0, 0, 0);
  }

  // Reading worked.
  return true;

}

//==========================================================================

// LHAupFromPYTHIA8 class.

//--------------------------------------------------------------------------

// Read in initialization information from PYTHIA 8.

bool LHAupFromPYTHIA8::setInit() {

  // Read in beam from Info class. Parton density left empty.
  int    idbmupA = infoPtr->idA();
  int    idbmupB = infoPtr->idB();
  double ebmupA  = infoPtr->eA();
  double ebmupB  = infoPtr->eB();
  int    pdfgupA = 0;
  int    pdfgupB = 0;
  int    pdfsupA = 0;
  int    pdfsupB = 0;
  setBeamA(idbmupA, ebmupA, pdfgupA, pdfsupA);
  setBeamB(idbmupB, ebmupB, pdfgupB, pdfsupB);

  // Currently only one allowed strategy.
  int    idwtup = 3;
  setStrategy(idwtup);

  // Only one process with dummy information. (Can overwrite at the end.)
  int    lprup  = 9999;
  double xsecup = 1.;
  double xerrup = 0.;
  double xmaxup = 1.;
  addProcess(lprup, xsecup, xerrup, xmaxup);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Read in event information from PYTHIA 8.

bool LHAupFromPYTHIA8::setEvent( int) {

  // Read process information from Info class, and store it.
  // Note: renormalization scale here, factorization further down.
  // For now always convert to process 9999, instead of infoPtr->code().
  int    idprup = 9999;
  double xwgtup = infoPtr->weight();
  double scalup = infoPtr->QRen();
  double aqedup = infoPtr->alphaEM();
  double aqcdup = infoPtr->alphaS();
  setProcess(idprup, xwgtup, scalup, aqedup, aqcdup);

  // Read in particle info one by one, excluding zero and beams, and store it.
  // Note unusual C++ loop range, to better reflect LHA/Fortran standard.
  int nup   = processPtr->size() - 3;
  int    idup, statusup, istup, mothup1, mothup2, icolup1, icolup2;
  double pup1, pup2, pup3, pup4, pup5, vtimup, spinup;
  for (int ip = 1; ip <= nup; ++ip) {
    Particle& particle = (*processPtr)[ip + 2];
    idup     = particle.id();
    // Convert from PYTHIA8 to LHA status codes.
    statusup = particle.status();
    if (ip < 3)            istup = -1;
    else if (statusup < 0) istup =  2;
    else                   istup =  1;
    mothup1  = max(0, particle.mother1() - 2);
    mothup2  = max(0, particle.mother2() - 2);
    icolup1  = particle.col();
    icolup2  = particle.acol();
    pup1     = particle.px();
    pup2     = particle.py();
    pup3     = particle.pz();
    pup4     = particle.e();
    pup5     = particle.m();
    vtimup   = particle.tau();
    spinup   = particle.pol();
    addParticle(idup, istup, mothup1, mothup2, icolup1, icolup2,
      pup1, pup2, pup3, pup4, pup5, vtimup, spinup, -1.) ;
  }

  // Extract hard-process initiator information from Info class, and store it.
  int    id1up      = infoPtr->id1();
  int    id2up      = infoPtr->id2();
  double x1up       = infoPtr->x1();
  double x2up       = infoPtr->x2();
  setIdX( id1up, id2up, x1up, x2up);

  // Also extract pdf information from Info class, and store it.
  int    id1pdfup   = infoPtr->id1pdf();
  int    id2pdfup   = infoPtr->id2pdf();
  double x1pdfup    = infoPtr->x1pdf();
  double x2pdfup    = infoPtr->x2pdf();
  double scalePDFup = infoPtr->QFac();
  double pdf1up     = infoPtr->pdf1();
  double pdf2up     = infoPtr->pdf2();
  setPdf(id1pdfup, id2pdfup, x1pdfup, x2pdfup, scalePDFup, pdf1up,
    pdf2up, true);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

//  Update cross-section information at the end of the run.

bool LHAupFromPYTHIA8::updateSigma() {

  // Read out information from PYTHIA 8 and send it in to LHA.
  double sigGen = CONVERTMB2PB * infoPtr->sigmaGen();
  double sigErr = CONVERTMB2PB * infoPtr->sigmaErr();
  setXSec(0, sigGen);
  setXErr(0, sigErr);

  // Done.
  return true;

}

//==========================================================================

} // end namespace Pythia8
