/*
 * Copyright (c) 2018 Jack Poulson <jack@hodgestar.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#ifndef QUOTIENT_TIMER_IMPL_H_
#define QUOTIENT_TIMER_IMPL_H_

#include <chrono>
#include <string>

#include "quotient/macros.hpp"

#include "quotient/timer.hpp"

namespace quotient {

inline Timer::Timer(const std::string& name) QUOTIENT_NOEXCEPT
: name_(name), running_(false), last_interval_seconds_(0), total_seconds_(0) { }

inline const std::string& Timer::Name() const QUOTIENT_NOEXCEPT {
    return name_;
}

inline void Timer::Start() QUOTIENT_NOEXCEPT {
  last_time_ = std::chrono::steady_clock::now();
  running_ = true;
}

inline double Timer::Stop() QUOTIENT_NOEXCEPT {
  QUOTIENT_ASSERT(running_, "The Timer was stopped when it was not running.");
  last_interval_seconds_ = SecondsSinceLastStart();
  total_seconds_ += last_interval_seconds_;
  running_ = false;
  return last_interval_seconds_;
}

inline double Timer::SecondsSinceLastStart() const QUOTIENT_NOEXCEPT {
  if (running_) {
    auto now = std::chrono::steady_clock::now();
    const std::chrono::duration<double> duration = now - last_time_;
    return duration.count();
  }
  return last_interval_seconds_;
}

inline double Timer::TotalSeconds() const QUOTIENT_NOEXCEPT {
  if (running_) {
    return total_seconds_ + SecondsSinceLastStart();
  }
  return total_seconds_;
}

inline void Timer::Reset(const std::string& name) QUOTIENT_NOEXCEPT {
  name_ = name;
  total_seconds_ = 0;
  last_interval_seconds_ = 0;
  running_ = false;
}

} // namespace quotient

#endif // ifndef QUOTIENT_TIMER_IMPL_H_
