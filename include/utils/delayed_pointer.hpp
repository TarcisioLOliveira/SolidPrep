/*
 *   Copyright (C) 2025 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of SolidPrep
 *
 *   SolidPrep is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   SolidPrep is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with SolidPrep.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#ifndef UTILS_DELAYED_POINTER_HPP
#define UTILS_DELAYED_POINTER_HPP

#include <algorithm>
#include <list>
#include <memory>
#include <utility>

namespace utils{

template<typename T>
class DelayedPointerView;

template<typename T>
class DelayedPointer{
    public:
    ~DelayedPointer(){
        for(auto p:this->views){
            *p = nullptr;
        }
    }
    DelayedPointer():ptr(nullptr){}

    template<typename... Args>
    DelayedPointer(Args... args):
        ptr(std::make_unique<T>(args...)){}
    DelayedPointer(std::unique_ptr<T> p):
        ptr(std::move(p)){}

    DelayedPointer(DelayedPointer&) = delete;
    DelayedPointer(DelayedPointer&& p):
        ptr(std::move(p.ptr)), views(std::move(p.views)){}
    DelayedPointer& operator=(DelayedPointer&) = delete;
    DelayedPointer& operator=(DelayedPointer&& p){
        this->ptr = std::move(p.ptr);
        this->views = std::move(p.views);
        return *this;
    }
    DelayedPointer& operator=(std::unique_ptr<T>&& p){
        this->ptr = std::move(p);
        for(auto p:this->views){
            *p = this->ptr.get();
        }
        return *this;
    }

    template<typename... Args>
    void construct(Args... args){
        this->ptr.reset(new T(args...));
        for(auto p:this->views){
            *p = this->ptr.get();
        }
    }
    void set(std::unique_ptr<T> p){
        this->ptr = std::move(p);
        for(auto p:this->views){
            *p = this->ptr.get();
        }
    }
    void destroy(){
        this->ptr.reset(nullptr);
        for(auto p:this->views){
            *p = nullptr;
        }
    }
    T* get(){
        return this->ptr.get();
    }
    T* get() const{
        return this->ptr.get();
    }

    DelayedPointerView<T> get_view(){
        return DelayedPointerView<T>(this);
    }
    DelayedPointerView<T> get_view() const{
        return DelayedPointerView<T>(this);
    }

    void register_view(DelayedPointerView<T>& vptr){
        this->views.push_back(&vptr.ptr);
    }
    void remove_view(DelayedPointerView<T>& vptr){
        auto addr = &vptr.ptr;
        auto pos(std::find(views.begin(), views.end(), addr));
        if(pos != this->views.end()){
            this->views.erase(pos);
        }
    }
    T& operator*(){
        return *this->ptr;
    }
    T& operator*() const{
        return *this->ptr;
    }
    T* operator->(){
        return this->ptr.get();
    }
    T* operator->() const{
        return this->ptr.get();
    }

    private:

    std::unique_ptr<T> ptr;
    std::list<T**> views;
};

template<typename T>
class DelayedPointerView{
    public:
    friend class DelayedPointer<T>;

    DelayedPointerView() = delete;
    DelayedPointerView(const DelayedPointerView& p):
        ptr(p.ptr), sptr(p.sptr){
        sptr->register_view(*this);
    }
    DelayedPointerView(DelayedPointerView&& p):
        ptr(p.ptr), sptr(p.sptr){
        sptr->register_view(*this);
    }
    DelayedPointerView(nullptr_t):
        ptr(nullptr), sptr(nullptr)
    {}
    DelayedPointerView(DelayedPointer<T>* sptr):
        ptr(sptr->get()), sptr(sptr)
    {
        sptr->register_view(*this);
    }
    ~DelayedPointerView(){
        if(this->sptr != nullptr){
            this->sptr->remove_view(*this);
        }
    }

    T* get(){
        return this->ptr;
    }
    T* get() const{
        return this->ptr;
    }

    T& operator*(){
        return *this->ptr;
    }
    T& operator*() const{
        return *this->ptr;
    }

    T* operator->(){
        return this->ptr;
    }
    T* operator->() const{
        return this->ptr;
    }

    DelayedPointerView& operator=(const DelayedPointerView& p){
        if(this->sptr != nullptr){
            this->sptr->remove_view(*this);
        }
        this->ptr = p.ptr;
        this->sptr = p.sptr;
        this->sptr->register_view(*this);
        return *this;
    }
    DelayedPointerView& operator=(DelayedPointerView&& p){
        if(this->sptr != nullptr){
            this->sptr->remove_view(*this);
        }
        this->ptr = p.ptr;
        this->sptr = p.sptr;
        this->sptr->register_view(*this);
        return *this;
    }

    private:
    T* ptr;
    DelayedPointer<T>* sptr;
};

}

#endif
