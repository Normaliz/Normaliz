//   Copyright (c)  2017  John Abbott,  Anna M. Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.

//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.

//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <https://www.gnu.org/licenses/>.


#include "CoCoA/SignalWatcher.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"

#include <iostream>

namespace CoCoA
{

  namespace // anonymous
  {

    volatile int/*sig_atomic_t*/ SignalReceived = 0; // GLOBAL VARIABLE!!!
    
  } // end of anonymous namespace


  SignalWatcher::SignalWatcher(int sig, void FnPtr(int)):
      mySig(sig)
  {
    myPrevSigactionPtr = nullptr;
    // struct sigaction sa;
    // if (FnPtr == nullptr)
    //   sa.sa_handler = &SetSignalReceived; // default CoCoA signal handler
    // else
    //   sa.sa_handler = FnPtr;
    // sigemptyset(&sa.sa_mask);
    // sa.sa_flags = SA_RESTART; // Restart functions if interrupted by handler

    // myPrevSigactionPtr = new struct sigaction;
    // if (sigaction(sig, &sa, myPrevSigactionPtr) != 0)
    // {
    //   delete myPrevSigactionPtr;
    //   CoCoA_THROW_ERROR("Unable to set signal handler", "SignalWatcher ctor");
    // }
  }


  void SignalWatcher::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "SignalWatcher(sig=" << mySig;
    if (IamActive())
      out << ')';
    out << ", DEACTIVATED)";
  }


  void SignalWatcher::myDeactivate() noexcept
  {
    if (!IamActive()) return;
    // sigaction(mySig, myPrevSigactionPtr, nullptr);
    // delete myPrevSigactionPtr;
    // myPrevSigactionPtr = nullptr;
  }


  SignalWatcher::~SignalWatcher()
  {
    myDeactivate();
  }


  std::ostream& operator<<(std::ostream& out, const SignalWatcher& SW)
  {
    if (!out) return out;  // short-cut for bad ostreams
    SW.myOutputSelf(out);
    return out;
  }


  int GetAndResetSignalReceived() noexcept
  {
    const int sig = SignalReceived;
    SignalReceived = 0;
    return sig;
  }

  void SetSignalReceived(int sig) noexcept
  {
    // SANITY CHECK on value of sig???
    CoCoA_ASSERT(sig != 0);
    SignalReceived = sig;
  }


  //-------------------------------------------------------
  // Must define this because of virtual mem fns.
  InterruptedBySignal::~InterruptedBySignal()
  {}

  void InterruptedBySignal::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "CoCoA::InterruptedBySignal(signal=" << mySignal
        << ", context=\"" << myContext << "\")";
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SignalWatcher.C,v 1.11 2021/10/07 12:38:13 abbott Exp $
// $Log: SignalWatcher.C,v $
// Revision 1.11  2021/10/07 12:38:13  abbott
// Summary: Added nissing include
//
// Revision 1.10  2021/10/07 12:13:09  abbott
// Summary: Added assertion to SetSignalReceived (must be non-zero)
//
// Revision 1.9  2021/02/10 19:40:00  abbott
// Summary: Added noexcept (sometimes instead of throw()) -- redmine 1572
//
// Revision 1.8  2020/06/17 15:49:27  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.7  2020/02/11 16:56:42  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.6  2020/02/11 16:12:19  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.5  2020/01/09 13:24:31  abbott
// Summary: Fixed redmine 1385 (memory leak)
//
// Revision 1.4  2019/03/19 11:07:07  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.3  2017/07/23 15:28:47  abbott
// Summary: Renamed SetInterruptSignalReceived to SetSignalReceived
//
// Revision 1.2  2017/07/22 13:01:02  abbott
// Summary: Added new exception class InterruptedBySignal; some cleaning
//
// Revision 1.1  2017/07/21 13:21:23  abbott
// Summary: Split olf interrupt into two ==> new file SignalWatcher; refactored interrupt and CpuTimeLimit
//
//
