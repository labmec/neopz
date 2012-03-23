//$Id: TPZGuiInterface.h,v 1.1 2010-03-30 14:37:06 fortiago Exp $

#ifndef TPZGuiInterfaceH
#define TPZGuiInterfaceH

#include "iostream"

/** This class implements a very simple interface from PZ kernel to GUI.
 * It implements for instance the capability of cancel the execution
 * The GUI must define a derived class which reimplements the messages and update methods
 * for better messages.
 * @author Tiago Forti
 * @since 2010, March 30
 */
class TPZGuiInterface{

	protected:

	/** Flag indicating if the thread is alive or if its canceling was requested
	 */
	bool fCanceled;

	/** Message attribute for UpdateGUI method */
	std::string fMessage;

	/** Progress bar attributes for UpdateGUI method */
	int fProgressBarPos, fProgressBarMaxPos, fProgressBarMinPos;

	public:

	/** Default constructor */
	TPZGuiInterface();

	/** Default destructor */
	virtual ~TPZGuiInterface();

	/** Updates the GUI with start messages
	 * This method must be reimplemented in derived classes for better messages
	 */
	virtual void Start();

	/** Updates the GUI
	 * This method must be reimplemented in derived classes for better messages
	 */
	virtual void UpdateCaption();

	/** Updates the GUI with ending messages
	 * This method must be reimplemented in derived classes for better messages
	 */
	virtual void End();

	/** Show an error message
	 * This method must be reimplemented in derived classes for better messages
	 */
	virtual void ShowErrorMessage(std::string message);

	/** Defines the process has been canceled */
	void SetKilled();

	/** Returns whether the process was canceled */
	bool AmIKilled();

	/** Message attribute for UpdateGUI method */
	std::string &Message(){
		return fMessage;
	}

	/** Progress bar attributes for UpdateGUI method */
	int &ProgressBarPos(){
		return fProgressBarPos;
	}

	/** Progress bar attributes for UpdateGUI method */
	int &ProgressBarMaxPos(){
		return fProgressBarMaxPos;
	}

    /** Progress bar attributes for UpdateGUI method */
	int &ProgressBarMinPos(){
		return fProgressBarMinPos;
	}

};

//---------------------------------------------------------------------------
#endif
