#ifndef CoCoA_SignalWatcher_H
#define CoCoA_SignalWatcher_H

//   Copyright (c)  2017  John Abbott,  Anna M. Bigatti

//   ALL RIGHTS RESERVED.  Open-source version: http://cocoa.dima.unige.it/


#include "CoCoA/interrupt.H"

#include <iosfwd>
//using std::ostream;
#include <csignal>
//using sigaction; (perhaps not in namespace std?)

namespace CoCoA
{

  class InterruptedBySignal: public InterruptReceived
  {
  public:
    InterruptedBySignal(int sig, const char* const context): InterruptReceived(context), mySignal(sig) {}
    virtual ~InterruptedBySignal();
    void myOutputSelf(std::ostream& out) const override;
  private: // data members (in addition to those of CoCoA::InterruptReceived)
    const int mySignal; // signal which caused the interrupt
    friend int TriggeringSignal(const InterruptedBySignal& intr) noexcept;
  };

  inline int TriggeringSignal(const InterruptedBySignal& intr) noexcept
  { return intr.mySignal; }


  // RAII class: ctor sets new sigaction, dtor restores previous one.
  // NB myDeactivate restores previous sigaction.
  class SignalWatcher
  {
  public:
    explicit SignalWatcher(int sig, void FnPtr(int) = nullptr);
    ~SignalWatcher();
    void myOutputSelf(std::ostream& out) const;
    friend void deactivate(SignalWatcher& SW) noexcept { SW.myDeactivate(); }
  private:
    bool IamActive() const noexcept { return (myPrevSigactionPtr != nullptr); }
    void myDeactivate() noexcept;
  private: // data members
    int mySig;
    void*   /* struct sigaction**/ myPrevSigactionPtr;
  };

  void deactivate(SignalWatcher& SW) noexcept;
  std::ostream& operator<<(std::ostream& out, const SignalWatcher& w);

  int GetAndResetSignalReceived() noexcept;
//  void ResetInterruptSignalReceived();
  void SetSignalReceived(int sig) noexcept; // used by interpreter & C5 (gui)

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SignalWatcher.H,v 1.7 2021/10/30 17:17:05 abbott Exp $
// $Log: SignalWatcher.H,v $
// Revision 1.7  2021/10/30 17:17:05  abbott
// Summary: Used keyword override (redmine 1625)
//
// Revision 1.6  2021/02/10 19:32:01  abbott
// Summary: Added noexcept (redmine 1572)
//
// Revision 1.5  2020/02/11 17:16:17  abbott
// Summary: All mem fns like myOutput now return void (redmine 969)
//
// Revision 1.4  2019/03/18 17:13:24  abbott
// Summary: Used nullptr instead of 0
//
// Revision 1.3  2017/07/23 15:28:46  abbott
// Summary: Renamed SetInterruptSignalReceived to SetSignalReceived
//
// Revision 1.2  2017/07/22 13:00:59  abbott
// Summary: Added new exception class InterruptedBySignal; some cleaning
//
// Revision 1.1  2017/07/21 13:21:22  abbott
// Summary: Split olf interrupt into two ==> new file SignalWatcher; refactored interrupt and CpuTimeLimit
//
//

#endif
