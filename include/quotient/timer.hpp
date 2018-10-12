/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_TIMER_H_
#define QUOTIENT_TIMER_H_

#include <chrono>
#include <string>

namespace quotient {

static constexpr char kDefaultTimerName[] = "[default]";

class Timer {
 public:
  // Constructs a timer with the given name.
  Timer(const std::string& name=kDefaultTimerName);

  // Returns the name of the timer.
  const std::string& Name() const;

  // (Re)starts the timer.
  void Start();

  // Stops the timer and returns the time in seconds since the last start.
  double Stop();

  // Returns the time in seconds since the last start.
  double SecondsSinceLastStart() const;

  // Returns the total time that the timer has run since the last 'Reset'.
  double TotalSeconds() const;

  // Resets the timer and changes its name to the specified value.
  void Reset(const std::string& name=kDefaultTimerName);

 private:
  // The name of the timer.
  std::string name_;

  // True if the timer is currently running.
  bool running_;

  // The time the timer was last started.
  std::chrono::steady_clock::time_point last_time_;

  // The number of seconds the timer has run since the last 'Start'.
  double last_interval_seconds_;

  // The total number of seconds the timer has run since the last 'Reset'.
  double total_seconds_; 
};

} // namespace quotient

#include "quotient/timer-impl.hpp"

#endif // ifndef QUOTIENT_TIMER_H_
