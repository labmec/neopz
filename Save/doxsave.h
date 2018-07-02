/**
 * @file
 * @brief Creates save group for Doxygen documentation.
 */

/**
 * @defgroup save Classes supporting persistency
 *
 * @brief Persistency within the PZ environment is implemented by deriving a class from the TPZSavable class and implementing the Read and Write method.
 *
 * Persistency (a.k.a. Object Serialization) in PZ has been improved in order to make backward object level compatibility possible. Thus, it
 * is possible to read objects of any class that were previously saved with its older versions. Class versioning is used to make sure old data
 * can be further used properly.
 *
 * For the improved serialization method to become effective, there are a few steps that need to be accomplished by a class developer when modifying a class attribute set
 * or creating a new class. 
 * \n 
 * \n Briefly, these are the necessary steps:
 * \n
 * \n 1) Write (or update if it already exists) a Translator class for a class that had changes in its attributes. This translator class must be able to read prior version objects saved in past.
 * \n 2) Write the method UpdateAttributes (if it doesn't exist) in Translator class.
 * \n 3) Write the method UpdateFromVx to make reading from last prior class version possible. This method will fit old data into the new class' attribute layout. 
 * \n   (x = last class version number prior to this new class version)
 * \n 4) Write (or update if it already exists) the method UpdateStream in Translator class. This method must call UpdateFromVx accordingly and UpdateAttributes.
 * \n
 */
