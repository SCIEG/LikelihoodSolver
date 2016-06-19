//------------------------------------------------------------------------------
// <auto-generated />
//
// This file was automatically generated by SWIG (http://www.swig.org).
// Version 3.0.10
//
// Do not make changes to this file unless you know what you are doing--modify
// the SWIG interface file instead.
//------------------------------------------------------------------------------


public class StringVectorVector : global::System.IDisposable, global::System.Collections.IEnumerable
    , global::System.Collections.Generic.IEnumerable<StringVector>
 {
  private global::System.Runtime.InteropServices.HandleRef swigCPtr;
  protected bool swigCMemOwn;

  internal StringVectorVector(global::System.IntPtr cPtr, bool cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = new global::System.Runtime.InteropServices.HandleRef(this, cPtr);
  }

  internal static global::System.Runtime.InteropServices.HandleRef getCPtr(StringVectorVector obj) {
    return (obj == null) ? new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero) : obj.swigCPtr;
  }

  ~StringVectorVector() {
    Dispose();
  }

  public virtual void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != global::System.IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          LabRetrieverPINVOKE.delete_StringVectorVector(swigCPtr);
        }
        swigCPtr = new global::System.Runtime.InteropServices.HandleRef(null, global::System.IntPtr.Zero);
      }
      global::System.GC.SuppressFinalize(this);
    }
  }

  public StringVectorVector(global::System.Collections.ICollection c) : this() {
    if (c == null)
      throw new global::System.ArgumentNullException("c");
    foreach (StringVector element in c) {
      this.Add(element);
    }
  }

  public bool IsFixedSize {
    get {
      return false;
    }
  }

  public bool IsReadOnly {
    get {
      return false;
    }
  }

  public StringVector this[int index]  {
    get {
      return getitem(index);
    }
    set {
      setitem(index, value);
    }
  }

  public int Capacity {
    get {
      return (int)capacity();
    }
    set {
      if (value < size())
        throw new global::System.ArgumentOutOfRangeException("Capacity");
      reserve((uint)value);
    }
  }

  public int Count {
    get {
      return (int)size();
    }
  }

  public bool IsSynchronized {
    get {
      return false;
    }
  }

  public void CopyTo(StringVector[] array)
  {
    CopyTo(0, array, 0, this.Count);
  }

  public void CopyTo(StringVector[] array, int arrayIndex)
  {
    CopyTo(0, array, arrayIndex, this.Count);
  }

  public void CopyTo(int index, StringVector[] array, int arrayIndex, int count)
  {
    if (array == null)
      throw new global::System.ArgumentNullException("array");
    if (index < 0)
      throw new global::System.ArgumentOutOfRangeException("index", "Value is less than zero");
    if (arrayIndex < 0)
      throw new global::System.ArgumentOutOfRangeException("arrayIndex", "Value is less than zero");
    if (count < 0)
      throw new global::System.ArgumentOutOfRangeException("count", "Value is less than zero");
    if (array.Rank > 1)
      throw new global::System.ArgumentException("Multi dimensional array.", "array");
    if (index+count > this.Count || arrayIndex+count > array.Length)
      throw new global::System.ArgumentException("Number of elements to copy is too large.");
    for (int i=0; i<count; i++)
      array.SetValue(getitemcopy(index+i), arrayIndex+i);
  }

  global::System.Collections.Generic.IEnumerator<StringVector> global::System.Collections.Generic.IEnumerable<StringVector>.GetEnumerator() {
    return new StringVectorVectorEnumerator(this);
  }

  global::System.Collections.IEnumerator global::System.Collections.IEnumerable.GetEnumerator() {
    return new StringVectorVectorEnumerator(this);
  }

  public StringVectorVectorEnumerator GetEnumerator() {
    return new StringVectorVectorEnumerator(this);
  }

  // Type-safe enumerator
  /// Note that the IEnumerator documentation requires an InvalidOperationException to be thrown
  /// whenever the collection is modified. This has been done for changes in the size of the
  /// collection but not when one of the elements of the collection is modified as it is a bit
  /// tricky to detect unmanaged code that modifies the collection under our feet.
  public sealed class StringVectorVectorEnumerator : global::System.Collections.IEnumerator
    , global::System.Collections.Generic.IEnumerator<StringVector>
  {
    private StringVectorVector collectionRef;
    private int currentIndex;
    private object currentObject;
    private int currentSize;

    public StringVectorVectorEnumerator(StringVectorVector collection) {
      collectionRef = collection;
      currentIndex = -1;
      currentObject = null;
      currentSize = collectionRef.Count;
    }

    // Type-safe iterator Current
    public StringVector Current {
      get {
        if (currentIndex == -1)
          throw new global::System.InvalidOperationException("Enumeration not started.");
        if (currentIndex > currentSize - 1)
          throw new global::System.InvalidOperationException("Enumeration finished.");
        if (currentObject == null)
          throw new global::System.InvalidOperationException("Collection modified.");
        return (StringVector)currentObject;
      }
    }

    // Type-unsafe IEnumerator.Current
    object global::System.Collections.IEnumerator.Current {
      get {
        return Current;
      }
    }

    public bool MoveNext() {
      int size = collectionRef.Count;
      bool moveOkay = (currentIndex+1 < size) && (size == currentSize);
      if (moveOkay) {
        currentIndex++;
        currentObject = collectionRef[currentIndex];
      } else {
        currentObject = null;
      }
      return moveOkay;
    }

    public void Reset() {
      currentIndex = -1;
      currentObject = null;
      if (collectionRef.Count != currentSize) {
        throw new global::System.InvalidOperationException("Collection modified.");
      }
    }

    public void Dispose() {
        currentIndex = -1;
        currentObject = null;
    }
  }

  public void Clear() {
    LabRetrieverPINVOKE.StringVectorVector_Clear(swigCPtr);
  }

  public void Add(StringVector x) {
    LabRetrieverPINVOKE.StringVectorVector_Add(swigCPtr, StringVector.getCPtr(x));
    if (LabRetrieverPINVOKE.SWIGPendingException.Pending) throw LabRetrieverPINVOKE.SWIGPendingException.Retrieve();
  }

  private uint size() {
    uint ret = LabRetrieverPINVOKE.StringVectorVector_size(swigCPtr);
    return ret;
  }

  private uint capacity() {
    uint ret = LabRetrieverPINVOKE.StringVectorVector_capacity(swigCPtr);
    return ret;
  }

  private void reserve(uint n) {
    LabRetrieverPINVOKE.StringVectorVector_reserve(swigCPtr, n);
  }

  public StringVectorVector() : this(LabRetrieverPINVOKE.new_StringVectorVector__SWIG_0(), true) {
  }

  public StringVectorVector(StringVectorVector other) : this(LabRetrieverPINVOKE.new_StringVectorVector__SWIG_1(StringVectorVector.getCPtr(other)), true) {
    if (LabRetrieverPINVOKE.SWIGPendingException.Pending) throw LabRetrieverPINVOKE.SWIGPendingException.Retrieve();
  }

  public StringVectorVector(int capacity) : this(LabRetrieverPINVOKE.new_StringVectorVector__SWIG_2(capacity), true) {
    if (LabRetrieverPINVOKE.SWIGPendingException.Pending) throw LabRetrieverPINVOKE.SWIGPendingException.Retrieve();
  }

  private StringVector getitemcopy(int index) {
    StringVector ret = new StringVector(LabRetrieverPINVOKE.StringVectorVector_getitemcopy(swigCPtr, index), true);
    if (LabRetrieverPINVOKE.SWIGPendingException.Pending) throw LabRetrieverPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  private StringVector getitem(int index) {
    StringVector ret = new StringVector(LabRetrieverPINVOKE.StringVectorVector_getitem(swigCPtr, index), false);
    if (LabRetrieverPINVOKE.SWIGPendingException.Pending) throw LabRetrieverPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  private void setitem(int index, StringVector val) {
    LabRetrieverPINVOKE.StringVectorVector_setitem(swigCPtr, index, StringVector.getCPtr(val));
    if (LabRetrieverPINVOKE.SWIGPendingException.Pending) throw LabRetrieverPINVOKE.SWIGPendingException.Retrieve();
  }

  public void AddRange(StringVectorVector values) {
    LabRetrieverPINVOKE.StringVectorVector_AddRange(swigCPtr, StringVectorVector.getCPtr(values));
    if (LabRetrieverPINVOKE.SWIGPendingException.Pending) throw LabRetrieverPINVOKE.SWIGPendingException.Retrieve();
  }

  public StringVectorVector GetRange(int index, int count) {
    global::System.IntPtr cPtr = LabRetrieverPINVOKE.StringVectorVector_GetRange(swigCPtr, index, count);
    StringVectorVector ret = (cPtr == global::System.IntPtr.Zero) ? null : new StringVectorVector(cPtr, true);
    if (LabRetrieverPINVOKE.SWIGPendingException.Pending) throw LabRetrieverPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public void Insert(int index, StringVector x) {
    LabRetrieverPINVOKE.StringVectorVector_Insert(swigCPtr, index, StringVector.getCPtr(x));
    if (LabRetrieverPINVOKE.SWIGPendingException.Pending) throw LabRetrieverPINVOKE.SWIGPendingException.Retrieve();
  }

  public void InsertRange(int index, StringVectorVector values) {
    LabRetrieverPINVOKE.StringVectorVector_InsertRange(swigCPtr, index, StringVectorVector.getCPtr(values));
    if (LabRetrieverPINVOKE.SWIGPendingException.Pending) throw LabRetrieverPINVOKE.SWIGPendingException.Retrieve();
  }

  public void RemoveAt(int index) {
    LabRetrieverPINVOKE.StringVectorVector_RemoveAt(swigCPtr, index);
    if (LabRetrieverPINVOKE.SWIGPendingException.Pending) throw LabRetrieverPINVOKE.SWIGPendingException.Retrieve();
  }

  public void RemoveRange(int index, int count) {
    LabRetrieverPINVOKE.StringVectorVector_RemoveRange(swigCPtr, index, count);
    if (LabRetrieverPINVOKE.SWIGPendingException.Pending) throw LabRetrieverPINVOKE.SWIGPendingException.Retrieve();
  }

  public static StringVectorVector Repeat(StringVector value, int count) {
    global::System.IntPtr cPtr = LabRetrieverPINVOKE.StringVectorVector_Repeat(StringVector.getCPtr(value), count);
    StringVectorVector ret = (cPtr == global::System.IntPtr.Zero) ? null : new StringVectorVector(cPtr, true);
    if (LabRetrieverPINVOKE.SWIGPendingException.Pending) throw LabRetrieverPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public void Reverse() {
    LabRetrieverPINVOKE.StringVectorVector_Reverse__SWIG_0(swigCPtr);
  }

  public void Reverse(int index, int count) {
    LabRetrieverPINVOKE.StringVectorVector_Reverse__SWIG_1(swigCPtr, index, count);
    if (LabRetrieverPINVOKE.SWIGPendingException.Pending) throw LabRetrieverPINVOKE.SWIGPendingException.Retrieve();
  }

  public void SetRange(int index, StringVectorVector values) {
    LabRetrieverPINVOKE.StringVectorVector_SetRange(swigCPtr, index, StringVectorVector.getCPtr(values));
    if (LabRetrieverPINVOKE.SWIGPendingException.Pending) throw LabRetrieverPINVOKE.SWIGPendingException.Retrieve();
  }

}
